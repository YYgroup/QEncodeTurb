"""Independent analysis of exported velocities; no import from geoencode_turb."""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import threading
import time
import numpy as np
import scipy
import psutil
from threadpoolctl import threadpool_limits,threadpool_info
from mps_ncs import to_scale_tensor,from_scale_tensor,spectra,rank_for_error,reconstruct
from mps_uniform import tt_svd

from paths import ROOT, CODE, NS, BLAS_THREADS
EPSILONS=(1e-2,1e-3,1e-4)

def stamp(**record):
    print(json.dumps(record),flush=True)

def worker(n,component):
    folder=ROOT/'mps_uniform'/f'n{n}'/component; folder.mkdir(parents=True,exist_ok=True)
    finished=folder/'complete.json'
    if finished.exists():
        saved=json.loads(finished.read_text())
        if not all((folder/a['file']).is_file() for a in saved['approximations']):
            raise RuntimeError('Checkpoint exists but an MPS core file is missing; use a fresh output directory.')
        input_hash=json.loads((ROOT/'fields'/f'n{n}'/'complete.json').read_text())['fields'][component]['sha256']
        if input_hash!=saved['input_sha256']:
            raise RuntimeError('Checkpoint input changed; use a fresh output directory.')
        if saved['adapter_sha256']!=hashlib.sha256((CODE/'mps_uniform.py').read_bytes()).hexdigest():
            raise RuntimeError('MPS adapter changed since completed calculation')
        stamp(event='reuse',n=n,component=component)
        return
    with threadpool_limits(limits=BLAS_THREADS):
        start=time.perf_counter()
        field=np.load(ROOT/'fields'/f'n{n}'/f'{component}.npy')
        tensor=to_scale_tensor(field)
        norm_squared=float(np.dot(tensor.ravel(),tensor.ravel()))
        # Validate index permutation independently of the SVD.
        inverse=from_scale_tensor(tensor)
        if not np.array_equal(field,inverse):
            raise RuntimeError('Scale encoding failed exact inverse check')
        del inverse,field; gc.collect()
        sstart=time.perf_counter()
        previous=ROOT/'mps'/f'n{n}'/component/'schmidt_spectra.npz'
        cached=previous if previous.exists() else ROOT/'survey'/f'n{n}_{component}.npz'
        if cached.exists():
            raw=np.load(cached)
            records=[]
            for cut in range(1,n//3):
                s=raw[f'cut{cut}']; p=s*s/norm_squared; positive=p[p>0]
                records.append(dict(cut=cut,values=s,entropy_bits=float(-np.dot(positive,np.log2(positive))),
                    norm_squared_from_spectrum=float(np.dot(s,s)),method=('Cached LAPACK eigvalsh of smaller Gram matrix' if min(8**cut,8**(n//3-cut))>1024 else 'Cached LAPACK singular values')))
            spectral_seconds=None
        else:
            records=spectra(tensor,lambda r:stamp(event='spectrum',n=n,component=component,cut=r['cut']))
            spectral_seconds=time.perf_counter()-sstart
        np.savez(folder/'schmidt_spectra.npz',**{f'cut{r["cut"]}':r['values'] for r in records})
        spectral_records=[]
        for r in records:
            spectral_records.append({key:value for key,value in r.items() if key!='values'}|
                dict(relative_norm_squared_residual=abs(r['norm_squared_from_spectrum']/norm_squared-1),
                    ranks_by_epsilon={str(e):rank_for_error(r['values'],e*e*norm_squared) for e in EPSILONS}))
        approximations=[]
        process=psutil.Process()
        for epsilon in EPSILONS:
            label=f'{epsilon:.0e}'
            path=folder/f'eps_{label}.json'
            if path.exists():
                approximations.append(json.loads(path.read_text())); continue
            baseline=process.memory_info().rss; peak=[baseline]; stop=threading.Event()
            def monitor():
                while not stop.wait(.01):
                    peak[0]=max(peak[0],process.memory_info().rss)
            thread=threading.Thread(target=monitor,daemon=True); thread.start()
            begin=time.perf_counter()
            cores,result=tt_svd(tensor,epsilon,[r['values'] for r in records],lambda cut,rank:stamp(event='tt_svd',n=n,component=component,epsilon=epsilon,cut=cut,rank=rank))
            build_seconds=time.perf_counter()-begin
            stop.set(); thread.join()
            rstart=time.perf_counter()
            recovered=reconstruct(cores)
            reconstruction_seconds=time.perf_counter()-rstart
            difference=(recovered-tensor).ravel()
            actual_error=float(np.sqrt(np.dot(difference,difference)/norm_squared))
            del difference
            if actual_error>epsilon*(1+1e-7)+1e-12:
                raise RuntimeError(f'Reconstruction exceeds error budget {actual_error}>{epsilon}')
            if abs(actual_error-result['certified_relative_error'])>1e-9:
                raise RuntimeError('Actual and certified TT-SVD errors disagree')
            physical=from_scale_tensor(recovered)
            middle=physical.shape[2]//2
            np.savez(folder/f'slices_eps_{label}.npz',
                reconstruction=physical[:,:,middle],reference=from_scale_tensor(tensor)[:,:,middle])
            del recovered,physical
            save_start=time.perf_counter()
            np.savez(folder/f'mps_eps_{label}.npz',**{f'core{i}':core for i,core in enumerate(cores)})
            save_seconds=time.perf_counter()-save_start
            lower_bound=max(r['ranks_by_epsilon'][str(epsilon)] for r in spectral_records)
            result.update(dict(epsilon=epsilon,actual_relative_l2_error=actual_error,
                chi_max=max(result['ranks']),chi_schmidt_lower_bound=lower_bound,
                chi_excess_over_lower_bound=max(result['ranks'])-lower_bound,
                dense_real_parameters=tensor.size,
                stored_parameter_ratio_to_dense=result['stored_real_parameters']/tensor.size,
                gauge_reduced_parameter_ratio_to_dense=result['gauge_reduced_real_parameters']/tensor.size,
                mps_core_bytes=sum(core.nbytes for core in cores),
                construction_seconds=build_seconds,reconstruction_seconds=reconstruction_seconds,
                save_seconds=save_seconds,construction_baseline_rss_bytes=baseline,
                construction_peak_rss_bytes=peak[0],construction_incremental_peak_rss_bytes=peak[0]-baseline,
                file=f'mps_eps_{label}.npz'))
            path.write_text(json.dumps(result,indent=2),encoding='utf-8')
            approximations.append(result)
            stamp(event='approximation_complete',n=n,component=component,epsilon=epsilon,chi=result['chi_max'],error=actual_error,seconds=build_seconds)
            del cores; gc.collect()
        metadata=dict(n=n,grid=2**(n//3),component=component,local_dimension=8,sites=n//3,
            encoding='Gourianov 2022 octal scale order, coarse to fine, omega=x_bit+2*y_bit+4*z_bit',
            adapter_sha256=hashlib.sha256((CODE/'mps_uniform.py').read_bytes()).hexdigest(),
            input_sha256=json.loads((ROOT/'fields'/f'n{n}'/'complete.json').read_text())['fields'][component]['sha256'],
            spectrum_seconds=spectral_seconds,spectrum_time_note='null denotes reused survey spectra',
            norm_squared=norm_squared,spectra=spectral_records,approximations=approximations,
            total_worker_seconds=time.perf_counter()-start,
            versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
            blas=threadpool_info())
        finished.write_text(json.dumps(metadata,indent=2),encoding='utf-8')
        stamp(event='component_complete',n=n,component=component)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--n',type=int); parser.add_argument('--component')
    args=parser.parse_args()
    if args.component:
        return worker(args.n,args.component)
    for n in ((args.n,) if args.n else NS):
        for component in ('ux','uy','uz'):
            folder=ROOT/'mps_uniform'/f'n{n}'/component; folder.mkdir(parents=True,exist_ok=True)
            stamp(event='launch',n=n,component=component)
            with (folder/'run.log').open('w',encoding='utf-8') as log:
                subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),'--n',str(n),'--component',component],
                    cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
            stamp(event='finished',n=n,component=component)

if __name__=='__main__':
    main()
