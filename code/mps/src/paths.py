"""Portable package paths and the selected grid list."""
import os
from pathlib import Path
CODE = Path(__file__).resolve().parent
PACKAGE = CODE.parent
PROJECT = PACKAGE / 'vendor'
ROOT = Path(os.environ.get('MPS_WORKDIR', str(PACKAGE / 'generated'))).resolve()
ROOT.mkdir(parents=True, exist_ok=True)
NS = tuple(int(n) for n in os.environ.get('MPS_N_VALUES', '15,18,21,24,27').split(','))
BLAS_THREADS = int(os.environ.get('MPS_BLAS_THREADS', '8'))
REFERENCE = PACKAGE / 'reference'
