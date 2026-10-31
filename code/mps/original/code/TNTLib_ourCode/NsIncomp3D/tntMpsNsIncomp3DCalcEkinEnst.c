/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2020
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 *
 * Calculates total kinetic energy and enstrophy of the 3D flow field
 *
 * \return No return value.
 */
void tntMpsNsIncomp3DCalcEkinEnst(double *Ekin, double *Enst, tntNetwork v1, tntNetwork v2, tntNetwork v3, tntNetwork d1, tntNetwork d2, tntNetwork d3, int chi, int orth_centre){
    
    /* Declare relevant variables */
    tntNetwork d1v2; tntNetwork d1v3;
    tntNetwork d2v1; tntNetwork d2v3;
    tntNetwork d3v1; tntNetwork d3v2;
    double err, N;
    int L;
    
    /* Initialise */
    err = 0;
    L = tntMpsLength(v1);
    N = pow(8,L);
    
    d1v2 = tntNetworkCopy(v2); d1v3 = tntNetworkCopy(v3);
    d2v1 = tntNetworkCopy(v1); d2v3 = tntNetworkCopy(v3);
    d3v1 = tntNetworkCopy(v1); d3v2 = tntNetworkCopy(v2);
    
    /* Calculate normalised kinetic energy */
    *Ekin = tntMpsSelfProduct(v1,orth_centre) + tntMpsSelfProduct(v2,orth_centre) + tntMpsSelfProduct(v3,orth_centre);
    *Ekin *= 0.5/N;
    
    /* Calculate wave-functions */
    err += tntMpsMpoProduct(d1v2,d1, chi);
    err += tntMpsMpoProduct(d1v3,d1, chi);
    err += tntMpsMpoProduct(d2v1,d2, chi);
    err += tntMpsMpoProduct(d2v3,d2, chi);
    err += tntMpsMpoProduct(d3v1,d3, chi);
    err += tntMpsMpoProduct(d3v2,d3, chi);
    
    /* Calculate normalised enstrophy */
    *Enst =  tntMpsSelfProduct(d2v3,-1) + tntMpsSelfProduct(d3v2,-1) - 2*(tntMpsMpsProduct(d2v3,d3v2)).re;
    *Enst += tntMpsSelfProduct(d3v1,-1) + tntMpsSelfProduct(d1v3,-1) - 2*(tntMpsMpsProduct(d3v1,d1v3)).re;
    *Enst += tntMpsSelfProduct(d1v2,-1) + tntMpsSelfProduct(d2v1,-1) - 2*(tntMpsMpsProduct(d1v2,d2v1)).re;
    *Enst *= 0.5/N;
    
    /* Now done. Free and exit (and ignore err – it's for debugging only anyway). */
    tntNetworkFree(&d1v2);
    tntNetworkFree(&d1v3);
    tntNetworkFree(&d2v1);
    tntNetworkFree(&d2v3);
    tntNetworkFree(&d3v1);
    tntNetworkFree(&d3v2);
    
}
