"""Cross-check the vectorized adapter against actual unmodified MATLAB code."""
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import scipy.linalg as la
from threadpoolctl import threadpool_limits
from mps_ncs import to_scale_tensor,reconstruct,from_scale_tensor,spectra

from paths import ROOT as root, BLAS_THREADS
original=loadmat(root/'validation/original_outputs.mat')
inputs=loadmat(root/'validation/original_inputs.mat')
field=inputs['F']; tensor=to_scale_tensor(field)
norm=np.linalg.norm(field)

def fixed_chi(tensor,chi):
    rank=1; remainder=tensor; cores=[]
    for site in range(tensor.ndim-1):
        u,s,vh=la.svd(remainder.reshape(rank*8,-1),full_matrices=False,check_finite=False)
        kept=min(chi,len(s))
        cores.append(u[:,:kept].reshape(rank,8,kept))
        remainder=s[:kept,None]*vh[:kept]
        rank=kept
    cores.append(remainder.reshape(rank,8,1))
    return from_scale_tensor(reconstruct(cores))

with threadpool_limits(limits=BLAS_THREADS):
    truncated=fixed_chi(tensor,16)
    records=spectra(tensor)
checks=dict(
    original_full_reconstruction_relative_error=float(np.linalg.norm(original['F_full']-field)/norm),
    original_truncated_relative_error=float(np.linalg.norm(original['F_truncated']-field)/norm),
    adapter_vs_original_truncated_relative_difference=float(np.linalg.norm(truncated-original['F_truncated'])/norm),
    marker_original_reconstruction_relative_error=float(np.linalg.norm(original['marker_reconstructed']-inputs['marker'])/np.linalg.norm(inputs['marker'])),
    maximum_relative_spectrum_difference=max(float(np.linalg.norm(r['values']-original['original_spectra'][r['cut']-1,:len(r['values'])])/norm) for r in records),
    maximum_entropy_difference_bits=max(abs(r['entropy_bits']-original['original_entropy'][0,r['cut']-1]) for r in records),
    original_memComp_full=float(original['mem_fraction_full'][0,0]),
    original_memComp_truncated=float(original['mem_fraction_truncated'][0,0]),
    original_seconds=float(original['original_seconds'][0,0]))
assert checks['original_full_reconstruction_relative_error']<1e-10
assert checks['adapter_vs_original_truncated_relative_difference']<1e-9
assert checks['maximum_relative_spectrum_difference']<1e-10
assert checks['maximum_entropy_difference_bits']<1e-9
assert checks['marker_original_reconstruction_relative_error']<1e-10
checks['passed']=True
(root/'validation/original_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
