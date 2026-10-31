/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2020
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*!
 * Calculates fidelity between fields [u1,u2,u3] and [v1,v2,v3] .
 */
double tntMpsNsIncomp3DCalcFidelity(tntNetwork u1, tntNetwork u2, tntNetwork u3,
                                  tntNetwork v1, tntNetwork v2, tntNetwork v3,
                                  int orth_centres){
    
    double normu, normv, fid;
    
    normu = tntMpsSelfProduct(u1,orth_centres) + tntMpsSelfProduct(u2,orth_centres) + tntMpsSelfProduct(u3,orth_centres);
    normv = tntMpsSelfProduct(v1,orth_centres) + tntMpsSelfProduct(v2,orth_centres) + tntMpsSelfProduct(v3,orth_centres);
    
    fid = tntComplexToReal(tntMpsMpsProduct(u1,v1)) + tntComplexToReal(tntMpsMpsProduct(u2,v2)) + tntComplexToReal(tntMpsMpsProduct(u3,v3));
    fid *= 1.0/sqrt(normu*normv);
    
    return fid;
    
}
