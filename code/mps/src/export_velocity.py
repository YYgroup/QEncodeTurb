"""Export one realization per grid, using the existing scientific workflow unchanged."""
import argparse
from dataclasses import asdict
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import threading
import time

from paths import ROOT, PROJECT, CODE, NS
sys.path.insert(0, str(PROJECT))

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(16*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def worker(n):
    import numpy as np
    import scipy
    import psutil
    import qiskit
    from qiskit import transpile
    from geoencode_turb.config import WorkflowConfig
    from geoencode_turb.quantum_encoding import run_sparse_quantum_encoding, build_3d_coupled_circuit
    from geoencode_turb.flow_reconstruction import reconstruct_flow_field
    folder = ROOT/'fields'/f'n{n}'
    folder.mkdir(parents=True, exist_ok=True)
    if (folder/'complete.json').exists():
        saved = json.loads((folder/'complete.json').read_text())
        if all(digest(folder/f'{key}.npy') == value['sha256'] for key, value in saved['fields'].items()):
            print('Verified existing fields', n, flush=True)
            return
        raise RuntimeError('Existing field hash mismatch')
    cfg = WorkflowConfig(nx=n//3, ny=n//3, nz=n//3,
                         seed_spin_up=2025, seed_spin_down=2026,
                         data_dir=folder/'compiler', figure_dir=folder/'figures')
    process = psutil.Process()
    baseline = process.memory_info().rss
    peak = [baseline]; stop = threading.Event()
    def monitor():
        while not stop.wait(.02):
            peak[0] = max(peak[0], process.memory_info().rss)
    thread = threading.Thread(target=monitor, daemon=True); thread.start()
    start = time.perf_counter()
    encoding = run_sparse_quantum_encoding(cfg, tail_tolerance=1e-12,
        locality_radius=3, prepare_plot_data=False, report_circuit_depth=False)
    encoding_seconds = time.perf_counter()-start
    params = encoding.amp_params
    diagnostics = asdict(encoding.compiler_diagnostics)
    circuits = []
    for seed in (2025, 2026):
        qc = build_3d_coupled_circuit(cfg.nx, cfg.ny, cfg.nz, params, seed=seed)
        tqc = transpile(qc, basis_gates=['u','cx'], optimization_level=3, seed_transpiler=2025)
        circuits.append(dict(seed=seed, logical_operations=dict(qc.count_ops()),
                             logical_depth=qc.depth(), basis_operations=dict(tqc.count_ops()),
                             basis_depth=tqc.depth()))
    active_amplitude = sum(int(abs(b)>0) + int(np.count_nonzero(np.abs(w)>1e-3)) for b,w in params)
    fitted_amplitude = sum(1+len(w) for b,w in params)
    phase_per_component = n*(n+1)//2
    reconstruction_start = time.perf_counter()
    flow = reconstruct_flow_field(encoding, cfg)
    reconstruction_seconds = time.perf_counter()-reconstruction_start
    del encoding; gc.collect()
    records = {}
    for name in ('ux','uy','uz'):
        field = np.ascontiguousarray(getattr(flow, name), dtype=np.float64)
        path = folder/f'{name}.npy'; np.save(path, field, allow_pickle=False)
        records[name] = dict(shape=list(field.shape), dtype=str(field.dtype),
            norm=float(np.linalg.norm(field.ravel())), mean=float(field.mean()),
            minimum=float(field.min()), maximum=float(field.max()),
            bytes=field.nbytes, sha256=digest(path))
        print('Exported', n, name, field.shape, flush=True)
        del field
    stop.set(); thread.join()
    metadata = dict(n_spatial_qubits=n, grid=cfg.grid_size, realization=1,
        config=asdict(cfg), tail_tolerance=1e-12, locality_radius=3,
        velocity_definition='Unmodified reconstruct_flow_field and compute_velocity(FDM)',
        axis_order='Existing flow.ux/uy/uz array axes, preserved without changes',
        encoding_seconds=encoding_seconds, reconstruction_seconds=reconstruction_seconds,
        export_total_seconds=time.perf_counter()-start,
        export_baseline_rss_bytes=baseline, export_peak_rss_bytes=peak[0],
        compiler=diagnostics, circuits=circuits,
        active_amplitude_parameters=active_amplitude,
        fitted_amplitude_slots=fitted_amplitude,
        phase_parameters_per_spin_component=phase_per_component,
        total_generative_scalar_parameters=active_amplitude+2*phase_per_component+1,
        generative_parameter_definition='Shared active amplitude angles + independent one-body and pairwise phase angles for both spin components + population ratio; counts arbitrary angles rather than RNG seeds. This circuit prepares a spinor, not velocity amplitudes. Classical decoding is required to obtain explicit velocity fields.',
        fields=records, versions=dict(python=platform.python_version(), numpy=np.__version__,
                                     scipy=scipy.__version__, qiskit=qiskit.__version__),
        scientific_source_sha256={p.relative_to(PROJECT).as_posix():digest(p) for p in
            [PROJECT/'geoencode_turb'/f for f in ('config.py','quantum_encoding.py','flow_reconstruction.py')]})
    (folder/'complete.json').write_text(json.dumps(metadata, indent=2, default=str), encoding='utf-8')
    print('Completed export', n, metadata['export_total_seconds'], flush=True)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--n', type=int)
    args=parser.parse_args()
    if args.n:
        worker(args.n)
    else:
        for n in NS:
            folder=ROOT/'fields'/f'n{n}'; folder.mkdir(parents=True, exist_ok=True)
            print('Starting', n, flush=True)
            with (folder/'export.log').open('w', encoding='utf-8') as log:
                subprocess.run([sys.executable,'-u',str(CODE/'export_velocity.py'),'--n',str(n)],
                    cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, check=True)
            print('Finished', n, flush=True)

if __name__=='__main__':
    main()
