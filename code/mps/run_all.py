"""Reproduce the velocity-field MPS study without paths outside this folder."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parent


def launch(script, arguments, root, env):
    logs = root/'logs'
    logs.mkdir(exist_ok=True)
    suffix = '_'.join(str(x).replace('--','') for x in arguments)
    log = logs/(Path(script).stem + ('_'+suffix if suffix else '')+'.log')
    command = [sys.executable, '-u', str(PACKAGE/'src'/script), *map(str, arguments)]
    print('Running:', script, *arguments, flush=True)
    start = time.perf_counter()
    with log.open('w', encoding='utf-8') as stream:
        result = subprocess.run(command, cwd=PACKAGE, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        print(log.read_text(encoding='utf-8', errors='replace')[-6000:], flush=True)
        raise RuntimeError(f'{script} failed. See {log}')
    print(f'Finished in {time.perf_counter()-start:.1f} s; log: {log}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', nargs='+', type=int, default=[15,18,21,24,27],
                        choices=[15,18,21,24,27], help='Spatial registers; default: the entire study.')
    parser.add_argument('--output', type=Path, default=PACKAGE/'generated')
    parser.add_argument('--threads', type=int, default=8, help='BLAS threads for the MPS analysis.')
    parser.add_argument('--tex', choices=['auto','required','off'], default='auto')
    parser.add_argument('--check', action='store_true', help='Fresh n=15,18 end-to-end numerical check, without the six-panel figure.')
    parser.add_argument('--plot-only', action='store_true', help='Redraw from the supplied reference plotting data; no field/MPS computation.')
    parser.add_argument('--verify-only', action='store_true', help='Check sources and supplied reference records; no heavy computation.')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    ns = [15,18] if args.check else sorted(set(args.n))
    if args.threads < 1:
        parser.error('--threads must be positive')
    env = os.environ.copy()
    env.update(MPS_WORKDIR=str(root), MPS_N_VALUES=','.join(map(str,ns)),
               MPS_BLAS_THREADS=str(args.threads), MPS_TEX=args.tex,
               MPLBACKEND='Agg', PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
    launch('verify_package.py', [], root, env)
    if args.verify_only:
        return
    if args.plot_only:
        shutil.copytree(PACKAGE/'reference/figure', root/'figure', dirs_exist_ok=True)
        shutil.copy2(PACKAGE/'templates/plot_mps_sixpanel.py', root/'figure/plot_mps_sixpanel.py')
        shutil.copy2(PACKAGE/'templates/plot_mps_sixpanel.ipynb', root/'figure/plot_mps_sixpanel.ipynb')
        env['MPS_N_VALUES'] = '15,18,21,24,27'
        launch('render_figure.py', [], root, env)
        print('Redrawn from supplied reference data:', root/'figure', flush=True)
        return

    # One launcher per output directory. A dead launcher's lock is recoverable.
    import psutil
    lock = root/'running.json'
    if lock.exists():
        previous = json.loads(lock.read_text(encoding='utf-8'))
        if psutil.pid_exists(previous['pid']):
            raise RuntimeError(f'Another launcher is using this output directory (PID {previous["pid"]}).')
        lock.unlink()
    with lock.open('x', encoding='utf-8') as stream:
        json.dump(dict(pid=os.getpid(), start=time.time()), stream)
    try:
        shutil.copytree(PACKAGE/'reference/validation', root/'validation', dirs_exist_ok=True)
        versions = {name: importlib.metadata.version(name) for name in
                    ['numpy','scipy','pandas','matplotlib','qiskit','qiskit-aer','psutil','threadpoolctl']}
        (root/'run_config.json').write_text(json.dumps(dict(n_values=ns,
            epsilons=[0.01,0.001,0.0001],threads=args.threads,versions=versions,
            source_manifest='mps/source_manifest.json',using_saved_matlab_outputs=True),indent=2),encoding='utf-8')
        launch('validate_original.py', [], root, env)
        for n in ns:
            launch('export_velocity.py', ['--n',n], root, env)
        launch('run_mps.py', [], root, env)
        launch('compile_circuit_30.py', [], root, env)
        launch('independent_recheck.py', [], root, env)
        for n in (15,27):
            if n in ns:
                launch('independent_recheck.py', ['--regenerate',n], root, env)
        launch('summarize.py', [], root, env)
        launch('prepare_plot_data.py', [], root, env)
        launch('compare_reference.py', [], root, env)
        if ns == [15,18,21,24,27]:
            launch('render_figure.py', [], root, env)
        else:
            print('Partial grids completed. The fixed six-panel publication figure requires all five grids.',flush=True)
        print('Study outputs:', root, flush=True)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
