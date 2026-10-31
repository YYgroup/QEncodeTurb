/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 *
 * Builds the advection-diffusion terms of the 2D incompressible Navier-Stokes equations
 *
 * \return No return value.
 */
void tntMpsNsIncomp2DBuildBurgerTerms(tntNetwork **BurgMpoSum, unsigned *numBurgTerms,
                                    tntNetwork u1, tntNetwork u2,
                                    tntNetwork d1, tntNetwork dd1, tntNetwork d2, tntNetwork dd2,
                                    double Re, double dt,int chi, int chi_mpo, double* errp){
    
    /* Declare the variables needed */
    /* In 2D, there are 6 Burger terms */
    *numBurgTerms = 6;
    *BurgMpoSum = malloc(6*sizeof(tntNetwork));
    
    /* Declare the terms that will make up BurgMpoSum */
    tntNetwork u1d1v;
    tntNetwork u2d2v;
    tntNetwork d1u1v;
    tntNetwork d2u2v;
    tntNetwork dd1v;
    tntNetwork dd2v;
        
    /* Declare dummy network */
    tntNetwork temp;
    
    /* Now fill up the terms */
    /* Start with advection terms */
    /* u1d1v */
    temp = tntMpsCreateHadamardMpo(u1);
    u1d1v = tntMpoMpoProduct(d1, temp);
    *errp += tntMpoTruncate(u1d1v, chi_mpo);
    
    /* d1u1v */
    d1u1v = tntMpoMpoProduct(temp, d1);
    *errp += tntMpoTruncate(d1u1v, chi_mpo);
    tntNetworkFree(&temp);
    
    /* u2d2v */
    temp = tntMpsCreateHadamardMpo(u2);
    u2d2v = tntMpoMpoProduct(d2, temp);
    *errp += tntMpoTruncate(u2d2v, chi_mpo);
    
    /* d2u2v */
    d2u2v = tntMpoMpoProduct(temp, d2);
    *errp += tntMpoTruncate(d2u2v, chi_mpo);
    tntNetworkFree(&temp);
    
    /* Now do the diffusion terms */
    /* dd1v, dd2v. */
    dd1v = tntNetworkCopy(dd1);
    dd2v = tntNetworkCopy(dd2);
    
    /* Scale the terms with relevant constants */
    tntNodeScaleReal(tntNodeFindFirst(u1d1v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(u2d2v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(d1u1v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(d2u2v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(dd1v),  +dt/Re);
    tntNodeScaleReal(tntNodeFindFirst(dd2v),  +dt/Re);
    
    /* Finally, put all the terms into the BurgMpoSum */
    /* BurgMpoSum: */
    (*BurgMpoSum)[0] = u1d1v;
    (*BurgMpoSum)[1] = u2d2v;
    (*BurgMpoSum)[2] = d1u1v;
    (*BurgMpoSum)[3] = d2u2v;
    (*BurgMpoSum)[4] = dd1v;
    (*BurgMpoSum)[5] = dd2v;
    
}
