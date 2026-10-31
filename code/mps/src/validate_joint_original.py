"""Validate the original joint-component MATLAB outputs independently."""
import json
import numpy as np
from scipy.io import loadmat
from scipy.linalg import svdvals
from threadpoolctl import threadpool_limits
from paths import ROOT, BLAS_THREADS
from mps_ncs import to_scale_tensor

def main():
    inputs=loadmat(ROOT/'recheck/joint_original_inputs.mat')
    output=loadmat(ROOT/'recheck/joint_original_outputs.mat')
    names=['U','V','W']; spectrum_names=['spectrumU','spectrumV','spectrumW']
    spectra=[]; norm_squared=0.; maximum_difference=0.
    with threadpool_limits(limits=BLAS_THREADS):
        for name,spectrum_name in zip(names,spectrum_names):
            field=inputs[name]; tensor=to_scale_tensor(field)
            norm=float(np.linalg.norm(field.ravel())); norm_squared+=norm**2
            cuts=[]
            for cut in range(1,tensor.ndim):
                values=svdvals(tensor.reshape(8**cut,-1),check_finite=False)
                original=output[spectrum_name][cut-1,:len(values)]
                maximum_difference=max(maximum_difference,float(np.linalg.norm(values-original)/norm))
                cuts.append(values)
            spectra.append(cuts)
    ranks=[]
    for epsilon in [0.01,0.001,0.0001]:
        per_cut=[]
        for cut in range(5):
            squared=sum(spectra[c][cut]**2 for c in range(3))
            tails=np.r_[np.cumsum(squared[::-1])[::-1],0.]
            per_cut.append(max(1,int(np.flatnonzero(tails<=epsilon**2*norm_squared)[0])))
        ranks.append(max(per_cut))
    assert ranks==np.max(output['chis'],axis=0).astype(int).tolist()
    assert maximum_difference<1e-10
    report=dict(passed=True,maximum_relative_spectrum_difference=maximum_difference,
                maximum_joint_cut_ranks=ranks,desired_relative_errors=[0.01,0.001,0.0001])
    (ROOT/'recheck/joint_original_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
