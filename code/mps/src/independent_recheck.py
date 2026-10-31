"""Independent numerical audit using physical coarse/fine blocks and saved cores.

This file deliberately does not import any MPS constructor, tensor permutation,
rank-selection function, or reconstructor from the previous calculation.
"""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
import scipy.linalg as la
from threadpoolctl import threadpool_limits

from paths import ROOT, PROJECT, NS, BLAS_THREADS
OUT=ROOT/'recheck'
EPS=(.01,.001,.0001)
COMPONENTS=('ux','uy','uz')

def emit(**values):
    print(json.dumps(values),flush=True)

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(16*1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def retained_rank(values,epsilon,norm_squared):
    # Independent forward prefix calculation in long double, rather than the
    # original reverse-cumulative-sum implementation.
    squared=np.asarray(values,dtype=np.longdouble)**2
    target=np.sum(squared,dtype=np.longdouble)-np.longdouble(epsilon)**2*np.longdouble(norm_squared)
    prefix=np.cumsum(squared,dtype=np.longdouble)
    return min(len(values),max(1,int(np.searchsorted(prefix,target,side='left'))+1))

def physical_unfold(field,cut):
    length=int(np.log2(field.shape[0])); coarse=2**cut; fine=2**(length-cut)
    # Three physical coordinates, each split into coarse and fine position.
    return field.reshape(coarse,fine,coarse,fine,coarse,fine).transpose(0,2,4,1,3,5).reshape(coarse**3,fine**3)

def contract_saved_cores(path,length):
    with np.load(path) as archive:
        cores=[archive[f'core{site}'] for site in range(length)]
    ranks=[1]
    for site,core in enumerate(cores):
        assert core.shape[0]==ranks[-1] and core.shape[1]==8
        ranks.append(core.shape[2])
    assert ranks[-1]==1
    prefix=cores[0].reshape(8,ranks[1])
    for site in range(1,length):
        prefix=(prefix@cores[site].reshape(ranks[site],8*ranks[site+1])).reshape(8**(site+1),ranks[site+1])
    return prefix[:,0],ranks

def physical_error(flat,field):
    # Map physical coordinates to base-8 scale indices directly, in small slabs.
    # This neither uses the original tensor permutation nor its inverse.
    size=field.shape[0]; length=int(np.log2(size))
    y=np.arange(size,dtype=np.int64)[None,:,None]
    z=np.arange(size,dtype=np.int64)[None,None,:]
    error_squared=np.longdouble(0); maximum=0.
    for start in range(0,size,4):
        stop=min(start+4,size)
        x=np.arange(start,stop,dtype=np.int64)[:,None,None]
        indices=np.zeros((stop-start,size,size),dtype=np.int64)
        for level in range(length):
            bit=length-1-level
            digit=((x>>bit)&1)+2*((y>>bit)&1)+4*((z>>bit)&1)
            indices=8*indices+digit
        difference=flat[indices]-field[start:stop]
        error_squared+=np.sum(difference.astype(np.longdouble)**2,dtype=np.longdouble)
        maximum=max(maximum,float(np.max(np.abs(difference))))
    return float(error_squared),maximum

def numerical_audit():
    OUT.mkdir(exist_ok=True)
    (OUT/'direct_spectra').mkdir(exist_ok=True)
    report=dict(cases=[],spectrum_checks=[],reconstruction_checks=[],cartesian_checks=[],joint_criterion=[])
    for n in NS:
        length=n//3; saved_spectra=[]; norms=[]
        exported=json.loads((ROOT/'fields'/f'n{n}'/'complete.json').read_text())
        for component in COMPONENTS:
            start=time.perf_counter()
            field_path=ROOT/'fields'/f'n{n}'/f'{component}.npy'
            assert sha256(field_path)==exported['fields'][component]['sha256']
            field=np.load(field_path)
            assert field.dtype==np.float64 and field.shape==(2**length,)*3 and np.isfinite(field).all()
            norm_squared=float(np.sum(field.astype(np.longdouble)**2,dtype=np.longdouble))
            original=json.loads((ROOT/'mps_uniform'/f'n{n}'/component/'complete.json').read_text())
            original_spectra=np.load(ROOT/'mps_uniform'/f'n{n}'/component/'schmidt_spectra.npz')
            direct=[]
            for cut in range(1,length):
                matrix=physical_unfold(field,cut)
                # Direct LAPACK SVD at every cut, including the large cuts
                # previously treated through Gram eigenvalues.
                s=la.svdvals(matrix,check_finite=False,overwrite_a=False)
                previous=original_spectra[f'cut{cut}']
                ranks={str(e):retained_rank(s,e,norm_squared) for e in EPS}
                prior=original['spectra'][cut-1]['ranks_by_epsilon']
                assert ranks==prior,(n,component,cut,ranks,prior)
                difference=float(np.linalg.norm(s-previous)/np.sqrt(norm_squared))
                assert difference<1e-7
                probabilities=s*s/norm_squared
                positive=probabilities[probabilities>0]
                entropy=float(-np.sum(positive*np.log2(positive)))
                entropy_difference=abs(entropy-original['spectra'][cut-1]['entropy_bits'])
                assert entropy_difference<1e-8
                report['spectrum_checks'].append(dict(n=n,component=component,cut=cut,
                    shape=list(matrix.shape),relative_spectrum_difference=difference,
                    entropy_difference_bits=entropy_difference,ranks=ranks,
                    norm_squared_residual=abs(float(np.dot(s,s))/norm_squared-1)))
                direct.append(s)
                del matrix
                emit(event='direct_svd_checked',n=n,component=component,cut=cut)
            np.savez(OUT/'direct_spectra'/f'n{n}_{component}.npz',**{f'cut{k+1}':s for k,s in enumerate(direct)})
            for epsilon in EPS:
                a=next(a for a in original['approximations'] if a['epsilon']==epsilon)
                flat,ranks=contract_saved_cores(ROOT/'mps_uniform'/f'n{n}'/component/a['file'],length)
                error_squared,max_error=physical_error(flat,field)
                error=float(np.sqrt(error_squared/norm_squared))
                assert error<=epsilon*(1+1e-7)+1e-12
                assert abs(error-a['actual_relative_l2_error'])<1e-10
                stored=8*sum(left*right for left,right in zip(ranks[:-1],ranks[1:]))
                assert stored==a['stored_real_parameters']
                assert stored-sum(r*r for r in ranks[1:-1])==a['gauge_reduced_real_parameters']
                report['reconstruction_checks'].append(dict(n=n,component=component,epsilon=epsilon,
                    independently_measured_relative_error=error,
                    difference_from_previous_error=abs(error-a['actual_relative_l2_error']),
                    maximum_absolute_pointwise_error=max_error,stored_parameters=stored))
                del flat; gc.collect()
                emit(event='saved_mps_checked',n=n,component=component,epsilon=epsilon,error=error)
            # A different, conventional three-site tensorization has local
            # dimension N, hence maximum rank N rather than the octal bound.
            cartesian=[]
            size=field.shape[0]
            for shape in ((size,size*size),(size*size,size)):
                s=la.svdvals(field.reshape(shape),check_finite=False)
                cartesian.append({str(e):retained_rank(s,e,norm_squared) for e in EPS})
            report['cartesian_checks'].append(dict(n=n,component=component,theoretical_maximum=size,
                                                  required_cut_ranks=cartesian))
            report['cases'].append(dict(n=n,component=component,field_hash_verified=True,
                elapsed_seconds=time.perf_counter()-start,norm_squared=norm_squared))
            saved_spectra.append(direct); norms.append(norm_squared)
            del field; gc.collect()
            (OUT/'independent_progress.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        for epsilon in EPS:
            cut_ranks=[]
            for cut in range(length-1):
                # The original findCritChi3D combines the three component
                # squared norms and uses desiredAcc=epsilon**2.
                squared=sum(saved_spectra[c][cut]**2 for c in range(3))
                cut_ranks.append(retained_rank(np.sqrt(squared),epsilon,sum(norms)))
            report['joint_criterion'].append(dict(n=n,epsilon=epsilon,cut_ranks=cut_ranks,
                                                 maximum_cut_rank=max(cut_ranks)))
    report['passed']=True
    report['summary']=dict(checked_fields=len(report['cases']),checked_spectra=len(report['spectrum_checks']),
        checked_saved_mps=len(report['reconstruction_checks']),
        maximum_relative_spectrum_difference=max(s['relative_spectrum_difference'] for s in report['spectrum_checks']),
        maximum_entropy_difference_bits=max(s['entropy_difference_bits'] for s in report['spectrum_checks']),
        maximum_reconstruction_error_difference=max(a['difference_from_previous_error'] for a in report['reconstruction_checks']))
    (OUT/'independent_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    emit(event='independent_audit_complete',**report['summary'])

def regenerate(n):
    sys.path.insert(0,str(PROJECT))
    from geoencode_turb.config import WorkflowConfig
    from geoencode_turb.quantum_encoding import run_sparse_quantum_encoding
    from geoencode_turb.flow_reconstruction import reconstruct_flow_field
    folder=OUT/'regenerated'/f'n{n}';folder.mkdir(parents=True,exist_ok=True)
    cfg=WorkflowConfig(nx=n//3,ny=n//3,nz=n//3,seed_spin_up=2025,seed_spin_down=2026,
        data_dir=folder/'compiler',figure_dir=folder/'figures')
    encoded=run_sparse_quantum_encoding(cfg,tail_tolerance=1e-12,locality_radius=3,
        prepare_plot_data=False,report_circuit_depth=False)
    flow=reconstruct_flow_field(encoded,cfg)
    result=dict(n=n,relative_errors={})
    for name in COMPONENTS:
        old=np.load(ROOT/'fields'/f'n{n}'/f'{name}.npy',mmap_mode='r')
        new=getattr(flow,name)
        difference=(new-old).ravel()
        error=float(np.linalg.norm(difference)/np.linalg.norm(old.ravel()))
        assert error<1e-10
        result['relative_errors'][name]=error
        del difference
    result['passed']=True
    (folder/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    emit(event='fresh_generation_checked',**result)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--regenerate',type=int)
    args=parser.parse_args()
    if args.regenerate:
        regenerate(args.regenerate)
    else:
        with threadpool_limits(limits=BLAS_THREADS):
            numerical_audit()
