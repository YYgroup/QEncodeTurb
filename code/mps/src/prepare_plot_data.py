"""Derive all figure inputs from fresh fields, MPS records and independent audit."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from paths import ROOT, NS, PACKAGE

EPS = (0.01,0.001,0.0001)

def main():
    destination=ROOT/'figure'
    data=destination/'data'
    data.mkdir(parents=True,exist_ok=True)
    audit=json.loads((ROOT/'recheck/independent_validation.json').read_text(encoding='utf-8'))
    assert audit['passed']
    summary=pd.read_csv(ROOT/'results/summary.csv')
    components=pd.read_csv(ROOT/'results/components.csv')
    circuits=pd.read_csv(ROOT/'results/circuit_resources.csv')
    spectra={}; metadata=[]; slices={}; velocities=[]; settings={}; bounds=[]; errors=[]
    for n in NS:
        length=n//3; maximum=8**(length//2)
        exported=json.loads((ROOT/'fields'/f'n{n}'/'complete.json').read_text(encoding='utf-8'))
        settings[str(n)]={k:v for k,v in exported['config'].items() if k not in ('data_dir','figure_dir')}
        settings[str(n)].update(tail_tolerance=1e-12,locality_radius=3)
        tails=np.zeros(maximum+1)
        for component in ('ux','uy','uz'):
            norm2=next(case['norm_squared'] for case in audit['cases'] if case['n']==n and case['component']==component)
            source=np.load(ROOT/'recheck/direct_spectra'/f'n{n}_{component}.npz')
            field=np.load(ROOT/'fields'/f'n{n}'/f'{component}.npy',mmap_mode='r')
            slices[f'n{n}_{component}']=np.array(field[:,:,field.shape[2]//2])
            velocities.append(dict(n=n,N=field.shape[0],component=component,rms=float(np.sqrt(norm2/field.size)),
                 slice_axis=2,slice_index=field.shape[2]//2,full_field_sha256=exported['fields'][component]['sha256']))
            del field
            for cut in range(1,length):
                values=source[f'cut{cut}']; key=f'n{n}_{component}_cut{cut}'
                spectra[key]=values
                check=next(x for x in audit['spectrum_checks'] if x['n']==n and x['component']==component and x['cut']==cut)
                p=values**2/norm2; positive=p[p>0]
                metadata.append(dict(n=n,N=2**length,L=length,component=component,cut=cut,key=key,
                    norm_squared=norm2,entropy_bits=float(-np.sum(positive*np.log2(positive))),
                    full_cut_dimension=min(8**cut,8**(length-cut)),
                    rank_eps_1e_02=check['ranks']['0.01'],rank_eps_1e_03=check['ranks']['0.001'],
                    rank_eps_1e_04=check['ranks']['0.0001']))
                squared=values**2
                tail=np.r_[np.cumsum(squared[::-1])[::-1],0.]
                tails[:len(tail)]=np.maximum(tails[:len(tail)],np.sqrt(tail/norm2))
        for epsilon in EPS:
            row=summary[(summary['n']==n)&np.isclose(summary['epsilon'],epsilon)].iloc[0]
            joint=next(x for x in audit['joint_criterion'] if x['n']==n and x['epsilon']==epsilon)
            bounds.append(dict(n=n,grid=2**length,epsilon=epsilon,octal_sites=length,octal_local_dimension=8,
                octal_maximum_chi=maximum,binary_maximum_chi=2**(n//2),cartesian_three_site_maximum_chi=2**length,
                componentwise_necessary_chi=int(row['chi_lower_bound_max']),verified_constructed_chi=int(row['chi_achieved_max']),
                necessary_fraction_of_maximum=row['chi_lower_bound_max']/maximum,
                constructed_fraction_of_maximum=row['chi_achieved_max']/maximum,
                original_joint_single_cut_chi=joint['maximum_cut_rank'],
                measured_maximum_component_error=row['maximum_component_relative_error']))
        errors.extend(dict(n=n,grid=2**length,chi=chi,fraction_of_maximum=chi/maximum,
                           necessary_maximum_component_error=tails[chi]) for chi in range(1,maximum+1))
    np.savez_compressed(data/'schmidt_spectra.npz',**spectra)
    np.savez_compressed(data/'velocity_slices.npz',**slices)
    pd.DataFrame(metadata).to_csv(data/'schmidt_metadata.csv',index=False)
    pd.DataFrame(velocities).to_csv(data/'velocity_metadata.csv',index=False)
    summary.to_csv(data/'mps_summary.csv',index=False)
    components.to_csv(data/'mps_components.csv',index=False)
    circuits.to_csv(data/'circuit_resources.csv',index=False)
    pd.DataFrame(bounds).to_csv(data/'upper_bounds_and_error_conventions.csv',index=False)
    pd.DataFrame(errors).to_csv(data/'unavoidable_error_curves.csv',index=False)
    # Also retain the two derived tables alongside the numerical audit.
    pd.DataFrame(bounds).to_csv(ROOT/'recheck/upper_bounds_and_error_conventions.csv',index=False)
    pd.DataFrame(errors).to_csv(ROOT/'recheck/unavoidable_error_curves.csv',index=False)
    (data/'field_settings.json').write_text(json.dumps(settings,indent=2),encoding='utf-8')
    (data/'audit_summary.json').write_text(json.dumps(dict(passed=True,independent_audit=audit['summary'],
         matlab_reference_comparison=json.loads((ROOT/'validation/original_validation.json').read_text(encoding='utf-8'))),indent=2),encoding='utf-8')
    record=json.loads((ROOT/'circuit_only/n30/complete.json').read_text(encoding='utf-8'))
    (data/'circuit_n30_compilation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    (data/'provenance.json').write_text(json.dumps(dict(paper_doi='10.1038/s43588-021-00181-1',
        original_code_archive_doi='10.25592/uhhfdm.14234',realization=1,phase_seeds=[2025,2026],
        n_values=list(NS),data_origin='Fresh outputs from this run; no reference numerical results copied.'),indent=2),encoding='utf-8')
    (data/'README.md').write_text('Fresh figure inputs from the generated velocity fields, saved MPS records and independent direct-SVD audit.\nThe velocity slices are retained for data completeness and are not plotted in the six-panel figure.\n',encoding='utf-8')
    manifest={p.relative_to(destination).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(data.iterdir()) if p.is_file()}
    (destination/'data_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(dict(grids=list(NS),spectra=len(metadata),approximations=len(components),
                          plotting_data_files=len(manifest)),indent=2))

if __name__=='__main__':
    main()
