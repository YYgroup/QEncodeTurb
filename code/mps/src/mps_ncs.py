"""Numerical adaptation of Gourianov's 2022 octal MPS encoding.

The original MATLAB sources are preserved in original/. This module replaces
the explicit binary counter by a tensor permutation and adds error-controlled
TT-SVD and independently computed, untruncated Schmidt spectra. The mathematical
encoding, coarse-to-fine ordering and local dimension 8 follow createMPS3D.m.
"""
import numpy as np
import scipy.linalg as la

def to_scale_tensor(field):
    length = int(np.log2(field.shape[0]))
    if field.shape != (2**length,)*3:
        raise ValueError('Expected a cubic power-of-two grid')
    # Original octal digit omega = x_bit + 2*y_bit + 4*z_bit.
    permutation = [axis*length+level for level in range(length) for axis in (2,1,0)]
    return np.ascontiguousarray(field.reshape((2,)*(3*length)).transpose(permutation)).reshape((8,)*length)

def from_scale_tensor(tensor):
    length=tensor.ndim
    permutation=[axis*length+level for level in range(length) for axis in (2,1,0)]
    return tensor.reshape((2,)*(3*length)).transpose(np.argsort(permutation)).reshape((2**length,)*3)

def rank_for_error(singular_values, absolute_error_squared):
    tails=np.r_[np.cumsum((singular_values**2)[::-1])[::-1],0.0]
    valid=np.flatnonzero(tails <= absolute_error_squared)
    return max(1,int(valid[0]))

def spectra(tensor, progress=None):
    norm_squared=float(np.dot(tensor.ravel(),tensor.ravel()))
    records=[]
    for cut in range(1,tensor.ndim):
        matrix=tensor.reshape(8**cut,-1)
        if min(matrix.shape)<=1024:
            values=la.svdvals(matrix,check_finite=False)
            method='LAPACK singular values'
        else:
            # Gram eigenvalues provide every squared Schmidt coefficient.
            gram=matrix@matrix.T if matrix.shape[0]<=matrix.shape[1] else matrix.T@matrix
            eigenvalues=la.eigvalsh(gram,check_finite=False,overwrite_a=True)[::-1]
            values=np.sqrt(np.maximum(eigenvalues,0))
            method='LAPACK eigvalsh of smaller Gram matrix'
        p=values**2/norm_squared
        positive=p[p>0]
        rec=dict(cut=cut,values=values,entropy_bits=float(-np.dot(positive,np.log2(positive))),
                 norm_squared_from_spectrum=float(np.dot(values,values)),method=method)
        records.append(rec)
        if progress:
            progress(rec)
    return records

def tt_svd(tensor,epsilon,progress=None):
    """Coarse-to-fine TT-SVD with an equal squared-error budget per bond."""
    length=tensor.ndim
    norm_squared=float(np.dot(tensor.ravel(),tensor.ravel()))
    budget=epsilon**2*norm_squared/(length-1)
    remainder=tensor.reshape(1,-1)
    rank=1; cores=[]; discarded=[]
    for site in range(length-1):
        matrix=remainder.reshape(rank*8,-1)
        u,s,vh=la.svd(matrix,full_matrices=False,check_finite=False,lapack_driver='gesdd')
        new_rank=rank_for_error(s,budget)
        cores.append(u[:,:new_rank].reshape(rank,8,new_rank))
        discarded.append(float(np.dot(s[new_rank:],s[new_rank:])))
        remainder=s[:new_rank,None]*vh[:new_rank]
        rank=new_rank
        del matrix,u,s,vh
        if progress:
            progress(site+1,rank)
    cores.append(remainder.reshape(rank,8,1))
    return cores,dict(discarded_squared_norms=discarded,
        certified_relative_error=float(np.sqrt(sum(discarded)/norm_squared)),
        ranks=[core.shape[2] for core in cores[:-1]],
        stored_real_parameters=sum(core.size for core in cores),
        gauge_reduced_real_parameters=sum(core.size for core in cores)-sum(core.shape[2]**2 for core in cores[:-1]))

def reconstruct(cores):
    result=cores[0][0]
    for core in cores[1:]:
        result=np.tensordot(result,core,axes=(-1,0))
    return result[...,0]
