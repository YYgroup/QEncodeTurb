"""Build auditable tables and static plots from the completed MPS records."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from paths import ROOT, NS
OUTPUT=ROOT/'results'; OUTPUT.mkdir(exist_ok=True)
FIGURES=ROOT/'figures'; FIGURES.mkdir(exist_ok=True)
COMPONENTS=('ux','uy','uz')
EPSILONS=(.01,.001,.0001)

def save_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def save_figure(fig,name):
    fig.savefig(FIGURES/f'{name}.png',dpi=220,bbox_inches='tight')
    fig.savefig(FIGURES/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)

def main():
    records={}; exports={}; details=[]; cuts=[]; rows=[]; resources=[]
    validations=[]
    for n in NS:
        exports[n]=json.loads((ROOT/'fields'/f'n{n}'/'complete.json').read_text())
        for component in COMPONENTS:
            metadata=json.loads((ROOT/'mps_uniform'/f'n{n}'/component/'complete.json').read_text())
            records[n,component]=metadata
            for cut in metadata['spectra']:
                for epsilon in EPSILONS:
                    cuts.append(dict(n=n,grid=metadata['grid'],component=component,cut=cut['cut'],
                        maximum_rank=min(8**cut['cut'],8**(n//3-cut['cut'])),
                        entropy_bits=cut['entropy_bits'],epsilon=epsilon,
                        schmidt_rank_lower_bound=cut['ranks_by_epsilon'][str(epsilon)],
                        spectrum_norm_relative_residual=cut['relative_norm_squared_residual']))
            for approximation in metadata['approximations']:
                details.append(dict(n=n,grid=metadata['grid'],component=component,**{
                    key:value for key,value in approximation.items()
                    if not isinstance(value,(dict,list))}))
                validations.append(approximation['actual_relative_l2_error']<=approximation['epsilon']*(1+1e-7)+1e-12)
                validations.append(abs(approximation['actual_relative_l2_error']-approximation['certified_relative_error'])<1e-9)
        export=exports[n]
        resources.append(dict(n=n,grid=export['grid'],
            generative_scalar_parameters=export['total_generative_scalar_parameters'],
            active_shared_amplitude_parameters=export['active_amplitude_parameters'],
            phase_parameters_both_components=2*export['phase_parameters_per_spin_component'],
            population_parameters=1,
            pre_transform_cx_gates_both_components=sum(c['basis_operations'].get('cx',0) for c in export['circuits']),
            pre_transform_u_gates_both_components=sum(c['basis_operations'].get('u',0) for c in export['circuits']),
            pre_transform_depth_max=max(c['basis_depth'] for c in export['circuits']),
            classical_compile_seconds=export['compiler']['compile_time_s'],
            classical_encoding_simulation_seconds=export['encoding_seconds'],
            classical_flow_reconstruction_seconds=export['reconstruction_seconds'],
            dense_velocity_bytes=sum(f['bytes'] for f in export['fields'].values()),
            endpoint_note='Circuit counts cover the existing pre-transform spinor circuits only. They exclude Gray decoding, Fourier transforms, nonlinear velocity evaluation, measurement, and evolution.'))
        for epsilon in EPSILONS:
            parts=[next(a for a in records[n,c]['approximations'] if a['epsilon']==epsilon) for c in COMPONENTS]
            norms=[records[n,c]['norm_squared'] for c in COMPONENTS]
            necessary_profiles=[]
            for component in COMPONENTS:
                necessary_profiles.append([1]+[s['ranks_by_epsilon'][str(epsilon)] for s in records[n,component]['spectra']]+[1])
            required_storage_lower_bound=sum(8*sum(a*b for a,b in zip(profile[:-1],profile[1:])) for profile in necessary_profiles)
            stored=sum(p['stored_real_parameters'] for p in parts)
            dense=3*2**n
            rows.append(dict(n=n,grid=2**(n//3),epsilon=epsilon,
                chi_lower_bound_max=max(p['chi_schmidt_lower_bound'] for p in parts),
                chi_achieved_max=max(p['chi_max'] for p in parts),
                chi_achieved_min=min(p['chi_max'] for p in parts),
                maximum_possible_chi=8**((n//3)//2),
                maximum_component_relative_error=max(p['actual_relative_l2_error'] for p in parts),
                vector_relative_l2_error=float(np.sqrt(sum(p['actual_relative_l2_error']**2*s for p,s in zip(parts,norms))/sum(norms))),
                stored_real_parameters=stored,
                gauge_reduced_real_parameters=sum(p['gauge_reduced_real_parameters'] for p in parts),
                dense_real_parameters=dense,
                stored_parameter_ratio_to_dense=stored/dense,
                gauge_reduced_parameter_ratio_to_dense=sum(p['gauge_reduced_real_parameters'] for p in parts)/dense,
                necessary_storage_parameter_lower_bound=required_storage_lower_bound,
                necessary_storage_ratio_lower_bound=required_storage_lower_bound/dense,
                mps_core_bytes=sum(p['mps_core_bytes'] for p in parts),
                mps_bytes_to_generative_description_bytes=stored/export['total_generative_scalar_parameters'],
                maximum_entropy_bits=max(s['entropy_bits'] for c in COMPONENTS for s in records[n,c]['spectra']),
                construction_seconds_sum=sum(p['construction_seconds'] for p in parts),
                reconstruction_seconds_sum=sum(p['reconstruction_seconds'] for p in parts),
                maximum_worker_peak_rss_bytes=max(p['construction_peak_rss_bytes'] for p in parts)))
    if not all(validations):
        raise RuntimeError('MPS reconstruction validation failed')
    save_csv(OUTPUT/'summary.csv',rows)
    save_csv(OUTPUT/'components.csv',details)
    save_csv(OUTPUT/'schmidt_cuts.csv',cuts)
    save_csv(OUTPUT/'circuit_resources.csv',resources)
    (OUTPUT/'summary.json').write_text(json.dumps(dict(summary=rows,circuit_resources=resources),indent=2),encoding='utf-8')
    extra=json.loads((ROOT/'circuit_only/n30/complete.json').read_text(encoding='utf-8'))['row']
    resources.append(extra)
    save_csv(OUTPUT/'circuit_resources.csv',resources)
    (OUTPUT/'summary.json').write_text(json.dumps(dict(summary=rows,circuit_resources=resources),indent=2),encoding='utf-8')
    audit=json.loads((ROOT/'recheck/independent_validation.json').read_text(encoding='utf-8'))
    validation=dict(passed=all(validations) and audit['passed'],completed_components=len(records),
                    completed_approximations=len(details),independent_audit=audit['summary'])
    (OUTPUT/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    print(json.dumps(validation,indent=2))

if __name__=='__main__':
    main()
