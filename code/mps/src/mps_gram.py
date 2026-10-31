"""The same fixed-cap TT-SVD, using orthogonal projections for large matrices.

No division by small singular values is needed. Large wide matrices use left
Gram eigenvectors. Large tall matrices use right eigenvectors followed by QR.
Full-rank splits use exact identity/QR factorizations. Every discarded norm is
measured from the actual projection residual, and final reconstruction is
checked independently by run_mps.py.
"""
import numpy as np
import scipy.linalg as la
from mps_uniform import choose_uniform_cap

def tt_svd(tensor,epsilon,exact_spectra,progress=None):
    norm_squared=float(np.dot(tensor.ravel(),tensor.ravel()))
    cap,bound=choose_uniform_cap(exact_spectra,epsilon,norm_squared)
    remainder=tensor.reshape(1,-1)
    rank=1; cores=[]; discarded=[]; methods=[]
    for site in range(tensor.ndim-1):
        matrix=remainder.reshape(rank*8,-1)
        m,n=matrix.shape; kept=min(cap,m,n)
        if min(m,n)<=1024:
            u,s,vh=la.svd(matrix,full_matrices=False,check_finite=False,lapack_driver='gesdd')
            left=u[:,:kept]
            right=s[:kept,None]*vh[:kept]
            lost=float(np.dot(s[kept:],s[kept:]))
            method='LAPACK SVD'
            del u,s,vh
        elif kept==min(m,n):
            if m<=n:
                left=np.eye(m)
                right=matrix
                method='Exact identity split (no truncation)'
            else:
                left,right=la.qr(matrix,mode='economic',check_finite=False)
                method='Exact QR split (no truncation)'
            lost=0.
        elif m<=n:
            gram=matrix@matrix.T
            values,vectors=la.eigh(gram,subset_by_index=(m-kept,m-1),check_finite=False,overwrite_a=True,driver='evr')
            left=np.ascontiguousarray(vectors[:,::-1])
            right=left.T@matrix
            residual=matrix-left@right
            lost=float(np.dot(residual.ravel(),residual.ravel()))
            method='Left Gram eigenspace and orthogonal projection'
            del gram,values,vectors,residual
        else:
            gram=matrix.T@matrix
            values,vectors=la.eigh(gram,subset_by_index=(n-kept,n-1),check_finite=False,overwrite_a=True,driver='evr')
            v=np.ascontiguousarray(vectors[:,::-1])
            projected=matrix@v
            left,r=la.qr(projected,mode='economic',check_finite=False)
            right=r@v.T
            residual=matrix-projected@v.T
            lost=float(np.dot(residual.ravel(),residual.ravel()))
            method='Right Gram eigenspace, projection and QR'
            del gram,values,vectors,v,projected,r,residual
        cores.append(left.reshape(rank,8,kept))
        discarded.append(lost); methods.append(method)
        remainder=right; rank=kept
        del matrix,left,right
        if progress:
            progress(site+1,rank)
    cores.append(remainder.reshape(rank,8,1))
    return cores,dict(uniform_cap=cap,schmidt_tail_sufficient_error_bound=bound,
        factorization_methods=methods,discarded_squared_norms=discarded,
        certified_relative_error=float(np.sqrt(sum(discarded)/norm_squared)),
        ranks=[core.shape[2] for core in cores[:-1]],
        stored_real_parameters=sum(core.size for core in cores),
        gauge_reduced_real_parameters=sum(core.size for core in cores)-sum(core.shape[2]**2 for core in cores[:-1]))
