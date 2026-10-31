"""Compile the 30-qubit row without simulating a statevector or velocity field."""
from dataclasses import asdict
import hashlib
import json
import platform
import sys
import time
import numpy as np
import qiskit
from qiskit import transpile
from paths import ROOT, PROJECT
sys.path.insert(0,str(PROJECT))
from geoencode_turb.config import WorkflowConfig
from geoencode_turb.quantum_encoding import compile_sparse_amplitude_parameters, build_3d_coupled_circuit

def main():
    folder = ROOT/'circuit_only/n30'
    folder.mkdir(parents=True,exist_ok=True)
    cfg = WorkflowConfig(nx=10,ny=10,nz=10,seed_spin_up=2025,seed_spin_down=2026,
                         data_dir=folder/'compiler',figure_dir=folder/'figures')
    begin = time.perf_counter()
    compiled = compile_sparse_amplitude_parameters(cfg,tail_tolerance=1e-12,locality_radius=3)
    params = compiled.amp_params
    circuits = []
    for seed in (2025,2026):
        circuit = build_3d_coupled_circuit(10,10,10,params,seed=seed)
        translated = transpile(circuit,basis_gates=['u','cx'],optimization_level=3,seed_transpiler=2025)
        circuits.append(dict(seed=seed,logical_operations=dict(circuit.count_ops()),
            logical_depth=circuit.depth(),basis_operations=dict(translated.count_ops()),
            basis_depth=translated.depth()))
    active = sum(int(abs(b)>0)+int(np.count_nonzero(np.abs(w)>1e-3)) for b,w in params)
    row = dict(n=30,grid=1024,generative_scalar_parameters=active+930+1,
        active_shared_amplitude_parameters=active,phase_parameters_both_components=930,
        population_parameters=1,
        pre_transform_cx_gates_both_components=sum(c['basis_operations'].get('cx',0) for c in circuits),
        pre_transform_u_gates_both_components=sum(c['basis_operations'].get('u',0) for c in circuits),
        pre_transform_depth_max=max(c['basis_depth'] for c in circuits),
        classical_compile_seconds=compiled.diagnostics.compile_time_s,
        classical_encoding_simulation_seconds=None,classical_flow_reconstruction_seconds=None,
        dense_velocity_bytes=None,
        endpoint_note='Circuit-only compilation; no statevector, velocity reconstruction, measurement or evolution.')
    result = dict(row=row,config=asdict(cfg),compiler=asdict(compiled.diagnostics),circuits=circuits,
        amplitude_parameters=[dict(bias=float(b),weights=np.asarray(w).tolist()) for b,w in params],
        tail_tolerance=1e-12,locality_radius=3,basis_gates=['u','cx'],optimization_level=3,
        seed_transpiler=2025,statevector_simulation_performed=False,flow_reconstruction_performed=False,
        elapsed_seconds=time.perf_counter()-begin,
        versions=dict(python=platform.python_version(),numpy=np.__version__,qiskit=qiskit.__version__),
        scientific_source_sha256={f'geoencode_turb/{name}':hashlib.sha256((PROJECT/'geoencode_turb'/name).read_bytes()).hexdigest()
                                  for name in ('config.py','quantum_encoding.py','flow_reconstruction.py')})
    (folder/'complete.json').write_text(json.dumps(result,indent=2,default=str),encoding='utf-8')
    print(json.dumps(row,indent=2))

if __name__ == '__main__':
    main()
