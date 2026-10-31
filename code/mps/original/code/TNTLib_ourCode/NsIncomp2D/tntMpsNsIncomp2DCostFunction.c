/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 *
 * Returns the value of the Burger's cost function. This is ||e_1||^2 +  ||e_2||^2  = ||(v1-(u1+ad1Sum))||^2 + ||(v2-(u2+ad2Sum))||^2
 *
 * \return No return value.
 */
double tntMpsNsIncomp2DBurgersCostFunction(tntNetwork v1, tntNetwork v2,
                                           tntNetwork u1, tntNetwork u2,
                                           tntNetwork* BurgMpoSum,
                                           unsigned numBurgTerms,
                                           int chi){
    
    double CF1 = 0;
    double CF2 = 0;
    double dummy;
    int L = tntMpsLength(v1);
    tntNetwork *Bu1, *Bu2;
    
    Bu1 = malloc(numBurgTerms*sizeof(tntNetwork));
    Bu2 = malloc(numBurgTerms*sizeof(tntNetwork));
    
    /* First contract v1^2+u1^2-2*v1*u1 and ditto for 2nd dimension */
    /* Start with the first dimension */
    CF1 += tntMpsSelfProduct(v1,-1);
    CF1 += tntMpsSelfProduct(u1,-1);
    CF1 += -2*tntComplexToReal(tntMpsMpsProduct(v1,u1));
    
    /* Now second dimension */
    CF2 += tntMpsSelfProduct(v2,-1);
    CF2 += tntMpsSelfProduct(u2,-1);
    CF2 += -2*tntComplexToReal(tntMpsMpsProduct(v2,u2));
    
    /* Contract the Burger terms with u1 and u2: */
    for(int i=0; i<numBurgTerms; i++){
        Bu1[i] = tntNetworkCopy(u1);
        Bu2[i] = tntNetworkCopy(u2);
        dummy = tntMpsMpoProduct(Bu1[i],BurgMpoSum[i],chi);
        dummy = tntMpsMpoProduct(Bu2[i],BurgMpoSum[i],chi);
    }
    
    /* Now contract the rest */
    /* Start with 2*u1*Bu1 - 2*v1*Bu1 and ditto for dimension 2 */
    for(int i=0; i<numBurgTerms; i++){
        /* Start with the first dimension */
        CF1 += -2*tntComplexToReal(tntMpsMpsProduct(v1,Bu1[i]));
        CF1 += 2*tntComplexToReal(tntMpsMpsProduct(u1,Bu1[i]));
        
        /* Now second dimension */
        CF2 += -2*tntComplexToReal(tntMpsMpsProduct(v2,Bu2[i]));
        CF2 += 2*tntComplexToReal(tntMpsMpsProduct(u2,Bu2[i]));
    }
    
    /* Now get Bu1*Bu1 and ditto for dimension 2 */
    for(int i=0; i<numBurgTerms; i++){
        for(int j=0; j<numBurgTerms; j++){
            /* Start with the first dimension */
            CF1 += tntComplexToReal(tntMpsMpsProduct(Bu1[i],Bu1[j]));
            
            /* Now second dimension */
            CF2 += tntComplexToReal(tntMpsMpsProduct(Bu2[i],Bu2[j]));
        }
    }
    
    /* Multiply by 4^-L to get the (quadratic) CF in units of per gridpoint. */
    CF1 *= pow((double) 1/4,(double) L);
    CF2 *= pow((double) 1/4,(double) L);
    
    /* Free Bu1 and Bu2 */
    for(int i=0; i<numBurgTerms; i++){
        tntNetworkFree(&(Bu1[i]));
        tntNetworkFree(&(Bu2[i]));
    }
    
    /* Done. Now add CF1 and CF2 together and return this result. */
    return CF1+CF2;
    
}

/*! \ingroup ns
 *
 * Returns the value of the divergence cost function. That is ||d1u1 + d2u2||^2 = d1u1d1u1 + d1u1d2u2 + d2u2d1u1 + d2u2d2u2
 multiplied by lambda.
 *
 * \return No return value.
 */

double tntMpsNsIncomp2DDivergenceCostFunction(tntNetwork v1, tntNetwork v2,
                                      tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d2xd2,
                                      int chi){
    
    double CF = 0;
    double error = 0;
    int L = tntMpsLength(v1);
    tntNetwork d1xd1v1, d2xd1v1, d2xd2v2, d1xd2v2;
    tntNetwork d2xd1 = d1xd2; //d1 and d2 are commuting operators
    
    /* First calculate d1v1d1v1 */
    d1xd1v1 = tntNetworkCopy(v1);
    error = tntMpsMpoProduct(d1xd1v1,d1xd1, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d1xd1v1,v1));
    
    /* Then calculate d2v2d1v1 */
    d2xd1v1 = tntNetworkCopy(v1);
    error = tntMpsMpoProduct(d2xd1v1,d2xd1, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d2xd1v1,v2));

    /* Then calculate d2v2d2v2 */
    d2xd2v2 = tntNetworkCopy(v2);
    error = tntMpsMpoProduct(d2xd2v2,d2xd2, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d2xd2v2,v2));
    
    /* Then calculate d1v1d2v2 */
    d1xd2v2 = tntNetworkCopy(v2);
    error = tntMpsMpoProduct(d1xd2v2,d1xd2, chi);
    CF += tntComplexToReal(tntMpsMpsProduct(d1xd2v2,v1));
    
    /* Finished. Now free all the generated networks and return the result. */
    tntNetworkFree(&d1xd1v1);
    tntNetworkFree(&d2xd1v1);
    tntNetworkFree(&d2xd2v2);
    tntNetworkFree(&d1xd2v2);
    
    return CF;
}
