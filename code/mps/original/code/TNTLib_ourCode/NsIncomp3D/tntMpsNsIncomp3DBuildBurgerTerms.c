/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 *
 * Builds the advection-diffusion terms of the 3D incompressible Navier-Stokes equations
 *
 * \return No return value.
 */
void tntMpsNsIncomp3DBuildBurgerTerms(tntNetwork **BurgMpoSum, unsigned *numBurgTerms,
                                    tntNetwork u1, tntNetwork u2, tntNetwork u3,
                                    tntNetwork d1, tntNetwork dd1,
                                    tntNetwork d2, tntNetwork dd2,
                                    tntNetwork d3, tntNetwork dd3,
                                    double Re, double dt,int chi, int chi_mpo, double* errp){
    
    /* Declare the variables needed */
    /* In 3D, there are 9 Burger terms */
    *numBurgTerms = 9;
    *BurgMpoSum = malloc(9*sizeof(tntNetwork));
    
    /* Declare the terms that will make up BurgMpoSum */
    tntNetwork u1d1v;
    tntNetwork u2d2v;
    tntNetwork u3d3v;
    tntNetwork d1u1v;
    tntNetwork d2u2v;
    tntNetwork d3u3v;
    tntNetwork dd1v;
    tntNetwork dd2v;
    tntNetwork dd3v;
    
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
    
    /* u3d3v */
    temp = tntMpsCreateHadamardMpo(u3);
    u3d3v = tntMpoMpoProduct(d3, temp);
    *errp += tntMpoTruncate(u3d3v, chi_mpo);
       
    /* d3u3v */
    d3u3v = tntMpoMpoProduct(temp, d3);
    *errp += tntMpoTruncate(d3u3v, chi_mpo);
    tntNetworkFree(&temp);
        
    /* Now do the diffusion terms */
    /* dd1v, dd2v, dd3v. */
    dd1v = tntNetworkCopy(dd1);
    dd2v = tntNetworkCopy(dd2);
    dd3v = tntNetworkCopy(dd3);

    /* Scale the terms with relevant constants */
    tntNodeScaleReal(tntNodeFindFirst(u1d1v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(u2d2v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(u3d3v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(d1u1v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(d2u2v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(d3u3v), -dt*0.5);
    tntNodeScaleReal(tntNodeFindFirst(dd1v), +1*dt/Re);
    tntNodeScaleReal(tntNodeFindFirst(dd2v), +1*dt/Re);
    tntNodeScaleReal(tntNodeFindFirst(dd3v), +1*dt/Re);
    
    /* Finally, put all the terms into the BurgMpoSum array */
    (*BurgMpoSum)[0] = u1d1v;
    (*BurgMpoSum)[1] = u2d2v;
    (*BurgMpoSum)[2] = u3d3v;
    (*BurgMpoSum)[3] = d1u1v;
    (*BurgMpoSum)[4] = d2u2v;
    (*BurgMpoSum)[5] = d3u3v;
    (*BurgMpoSum)[6] = dd1v;
    (*BurgMpoSum)[7] = dd2v;
    (*BurgMpoSum)[8] = dd3v;
}
