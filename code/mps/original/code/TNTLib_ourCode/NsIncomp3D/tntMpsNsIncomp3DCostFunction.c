/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*!
 * Returns the value of the Burger's cost function. This is ||e_1||^2 +  ||e_2||^2 +  ||e_3||^2  = ||(v1-(u1+Bu1))||^2 + ||(v2-(u2+Bu2))||^2 +  ||(v3-(u3+Bu3))||^2
 */
double tntMpsNsIncomp3DBurgersCostFunction(tntNetwork v1, tntNetwork v2, tntNetwork v3,
                                   tntNetwork u1, tntNetwork u2, tntNetwork u3,
                                   tntNetwork* BurgMpoSum,
                                   unsigned numBurgTerms){
    
    double CF1 = 0;
    double CF2 = 0;
    double CF3 = 0;
    double dummy;
    int L = tntMpsLength(v1);
    tntNetwork *Bu1, *Bu2, *Bu3;
    
    Bu1 = malloc(numBurgTerms*sizeof(tntNetwork));
    Bu2 = malloc(numBurgTerms*sizeof(tntNetwork));
    Bu3 = malloc(numBurgTerms*sizeof(tntNetwork));
    
    /* First contract v1^2+u1^2-2*v1*u1 and ditto for 2nd and 3rd dimensions */
    /* Start with the first dimension */
    CF1 += tntMpsSelfProduct(v1,-1);
    CF1 += tntMpsSelfProduct(u1,-1);
    CF1 += -2*tntComplexToReal(tntMpsMpsProduct(v1,u1));
    
    /* Now second dimension */
    CF2 += tntMpsSelfProduct(v2,-1);
    CF2 += tntMpsSelfProduct(u2,-1);
    CF2 += -2*tntComplexToReal(tntMpsMpsProduct(v2,u2));
    
    /* Third dimension */
    CF3 += tntMpsSelfProduct(v3,-1);
    CF3 += tntMpsSelfProduct(u3,-1);
    CF3 += -2*tntComplexToReal(tntMpsMpsProduct(v3,u3));
    
    /* Contract the Burger terms with u1, u2 and u3: */
    for(int i=0; i<numBurgTerms; i++){
        Bu1[i] = tntNetworkCopy(u1);
        Bu2[i] = tntNetworkCopy(u2);
        Bu3[i] = tntNetworkCopy(u3);
        dummy = tntMpsMpoProduct(Bu1[i],BurgMpoSum[i],-1);
        dummy = tntMpsMpoProduct(Bu2[i],BurgMpoSum[i],-1);
        dummy = tntMpsMpoProduct(Bu3[i],BurgMpoSum[i],-1);
    }
    
    /* Now contract the rest */
    /* Start with 2*u1*Bu1 - 2*v1*Bu1 and ditto for dimensions 2 and 3 */
    for(int i=0; i<numBurgTerms; i++){
        /* Start with the first dimension */
        CF1 += -2*tntComplexToReal(tntMpsMpsProduct(v1,Bu1[i]));
        CF1 += 2*tntComplexToReal(tntMpsMpsProduct(u1,Bu1[i]));
        
        /* Now second dimension */
        CF2 += -2*tntComplexToReal(tntMpsMpsProduct(v2,Bu2[i]));
        CF2 += 2*tntComplexToReal(tntMpsMpsProduct(u2,Bu2[i]));
        
        /* Finally third dimension */
        CF3 += -2*tntComplexToReal(tntMpsMpsProduct(v3,Bu3[i]));
        CF3 += 2*tntComplexToReal(tntMpsMpsProduct(u3,Bu3[i]));
    }
    
    /* Now get Bu1*Bu1 and ditto for dimensions 2 and 3 */
    for(int i=0; i<numBurgTerms; i++){
        for(int j=0; j<numBurgTerms; j++){
            /* Start with the first dimension */
            CF1 += tntComplexToReal(tntMpsMpsProduct(Bu1[i],Bu1[j]));
            
            /* Now second dimension */
            CF2 += tntComplexToReal(tntMpsMpsProduct(Bu2[i],Bu2[j]));
            
            /* Now third dimension */
            CF3 += tntComplexToReal(tntMpsMpsProduct(Bu3[i],Bu3[j]));
        }
    }
    
    /* Free Bu1, Bu2 and Bu3 */
    for(int i=0; i<numBurgTerms; i++){
        tntNetworkFree(&(Bu1[i]));
        tntNetworkFree(&(Bu2[i]));
        tntNetworkFree(&(Bu3[i]));
    }
    
    /* Done. Now add CF1, CF2 and CF3 together and return this result. */
    return CF1+CF2+CF3;
    
}

/*!
* Returns the value of the divergence cost function. That is ||d1u1 + d2u2  + d3u3||^2 =  <u1d1d1u1> + <u2d2d1u1> + .... + <u3d3d3u3>
*/

double tntMpsNsIncomp3DDivergenceCostFunction(tntNetwork v1, tntNetwork v2, tntNetwork v3,
                                      tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d1xd3, tntNetwork d2xd2, tntNetwork d2xd3, tntNetwork d3xd3, int chi){
    
    double CF = 0;
    double error = 0;
    tntNetwork d1xd1v1, d2xd1v1, d3xd1v1, d1xd2v2, d2xd2v2, d3xd2v2, d1xd3v3, d2xd3v3, d3xd3v3;
    tntNetwork d2xd1, d3xd1, d3xd2;
    d2xd1 = d1xd2; d3xd1 = d1xd3; d3xd2 = d2xd3; //d1, d2 and d3 are commuting operators
    
    /* First calculate d1v1d1v1 */
    d1xd1v1 = tntNetworkCopy(v1);
    error = tntMpsMpoProduct(d1xd1v1,d1xd1, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d1xd1v1,v1));

    /* Then calculate <v2d2xd1v1> */
    d2xd1v1 = tntNetworkCopy(v1);
    error = tntMpsMpoProduct(d2xd1v1,d2xd1, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d2xd1v1,v2));

    /* Then calculate <v1d3xd1v1> */
    d3xd1v1 = tntNetworkCopy(v1);
    error = tntMpsMpoProduct(d3xd1v1,d3xd1, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d3xd1v1,v3));

    /* Now calculate <v1d1xd2v2> */
    d1xd2v2 = tntNetworkCopy(v2);
    error = tntMpsMpoProduct(d1xd2v2,d1xd2, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d1xd2v2,v1));

    /* Then calculate <v2d2xd2v2> */
    d2xd2v2 = tntNetworkCopy(v2);
    error = tntMpsMpoProduct(d2xd2v2,d2xd2, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d2xd2v2,v2));

    /* Then calculate <v3d3xd2v2> */
    d3xd2v2 = tntNetworkCopy(v2);
    error = tntMpsMpoProduct(d3xd2v2,d3xd2, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d3xd2v2,v3));

    /* Now calculate <v1d1xd3v3> */
    d1xd3v3 = tntNetworkCopy(v3);
    error = tntMpsMpoProduct(d1xd3v3,d1xd3, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d1xd3v3,v1));

    /* Then calculate <v2d2xd3v3> */
    d2xd3v3 = tntNetworkCopy(v3);
    error = tntMpsMpoProduct(d2xd3v3,d2xd3, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d2xd3v3,v2));

    /* Then calculate <v3d3xd3v3> */
    d3xd3v3 = tntNetworkCopy(v3);
    error = tntMpsMpoProduct(d3xd3v3,d3xd3, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d3xd3v3,v3));
    
   // /* Multiply by 8^-L to get the (quadratic) CF in units of per gridpoint. */
   // CF *= h*h*h*h;
    
    /* Finished. Now free all the generated networks and return the result. */
    tntNetworkFree(&d1xd1v1);
    tntNetworkFree(&d2xd1v1);
    tntNetworkFree(&d3xd1v1);
    tntNetworkFree(&d1xd2v2);
    tntNetworkFree(&d2xd2v2);
    tntNetworkFree(&d3xd2v2);
    tntNetworkFree(&d1xd3v3);
    tntNetworkFree(&d2xd3v3);
    tntNetworkFree(&d3xd3v3);
    
    return CF;
}
