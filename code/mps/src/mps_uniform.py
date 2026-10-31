"""Fixed-cap TT-SVD following the original NCS MPS constructor.

The cap is selected from the sum of the exact unfolding tail energies, which
provides a sufficient global error bound for TT-SVD. This removes the excess
rank from an equal local error budget when only one or two cuts dominate.
"""
import numpy as np
import scipy.linalg as la

def choose_uniform_cap(spectra,epsilon,norm_squared):
    maximum=max(len(s) for s in spectra)
    tail_sum=np.zeros(maximum+1)
    for s in spectra:
        tails=np.r_[np.cumsum((s*s)[::-1])[::-1],0.]
        tail_sum[:len(tails)]+=tails
    eligible=np.flatnonzero(tail_sum <= epsilon*epsilon*norm_squared)
    cap=max(1,int(eligible[0]))
    return cap,float(np.sqrt(tail_sum[cap]/norm_squared))

def tt_svd(tensor,epsilon,exact_spectra,progress=None):
    norm_squared=float(np.dot(tensor.ravel(),tensor.ravel()))
    cap,bound=choose_uniform_cap(exact_spectra,epsilon,norm_squared)
    remainder=tensor.reshape(1,-1)
    rank=1; cores=[]; discarded=[]
    for site in range(tensor.ndim-1):
        matrix=remainder.reshape(rank*8,-1)
        u,s,vh=la.svd(matrix,full_matrices=False,check_finite=False,lapack_driver='gesdd')
        kept=min(cap,len(s))
        cores.append(u[:,:kept].reshape(rank,8,kept))
        discarded.append(float(np.dot(s[kept:],s[kept:])))
        remainder=s[:kept,None]*vh[:kept]
        rank=kept
        del matrix,u,s,vh
        if progress:
            progress(site+1,rank)
    cores.append(remainder.reshape(rank,8,1))
    return cores,dict(uniform_cap=cap,schmidt_tail_sufficient_error_bound=bound,
        discarded_squared_norms=discarded,
        certified_relative_error=float(np.sqrt(sum(discarded)/norm_squared)),
        ranks=[core.shape[2] for core in cores[:-1]],
        stored_real_parameters=sum(core.size for core in cores),
        gauge_reduced_real_parameters=sum(core.size for core in cores)-sum(core.shape[2]**2 for core in cores[:-1]))
