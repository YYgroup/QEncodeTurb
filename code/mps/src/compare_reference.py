"""Compare scientific quantities with supplied records, excluding timings/RSS."""
import json
import numpy as np
import pandas as pd
from paths import ROOT, REFERENCE, NS

def main():
    reference=pd.read_csv(REFERENCE/'results/summary.csv')
    fresh=pd.read_csv(ROOT/'results/summary.csv')
    exact=['chi_lower_bound_max','chi_achieved_max','maximum_possible_chi',
           'stored_real_parameters','gauge_reduced_real_parameters','dense_real_parameters','mps_core_bytes']
    checks=[]
    for _, row in fresh.iterrows():
        prior=reference[(reference['n']==row['n'])&np.isclose(reference['epsilon'],row['epsilon'])].iloc[0]
        for key in exact:
            assert int(row[key])==int(prior[key]), (row['n'],row['epsilon'],key,row[key],prior[key])
        assert abs(row['maximum_component_relative_error']-prior['maximum_component_relative_error'])<1e-8
        checks.append(dict(n=int(row['n']),epsilon=float(row['epsilon']),passed=True))
    # Circuit counts depend on the pinned Qiskit version; make changes visible.
    circuits=pd.read_csv(ROOT/'results/circuit_resources.csv')
    prior=pd.read_csv(REFERENCE/'figure/data/circuit_resources.csv')
    for _, row in circuits.iterrows():
        r=prior[prior['n']==row['n']].iloc[0]
        for key in ['generative_scalar_parameters','pre_transform_cx_gates_both_components',
                    'pre_transform_u_gates_both_components','pre_transform_depth_max']:
            assert int(row[key])==int(r[key]), (row['n'],key,row[key],r[key])
    audit=json.loads((ROOT/'recheck/independent_validation.json').read_text(encoding='utf-8'))
    assert audit['passed']
    result=dict(passed=True,n_values=list(NS),checked_summary_rows=len(checks),
                checked_circuit_rows=len(circuits),tolerance_for_reconstruction_error_difference=1e-8,
                timing_and_memory_not_required_to_match=True,checks=checks)
    (ROOT/'results/reference_comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
