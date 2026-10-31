/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2020
*/

/* Include header for NS and PDE functions, as well as headers for some internal functions of the TNT library */
#include "tntPde.h"
#include "tntIncompNs.h"
#include "../lib/headers/tnt_int.h"
#include "../lib/headers/dec_pde_int.h"
#include "../lib/headers/dec_matrix.h"
#include "../lib/headers/dec_la.h"

/*! \ingroup ns
 * This function, along with the two others in this .c file, builds a preconditioner for use in the conjugate gradient descent optimisation algorithm of tntNodeNsIncomp2DConjGrad(...). Note: I have personally not  observed any significant speedup when using this preconditioner. Indeed, in my experience, it just slows the code down and therefore I reccommend against activating this preconditioner in tntNodeNsIncomp2DConjGrad(...).
 */

tntNode tntNodeNsIncomp2DConjGrad_BuildLeftEnvPreconditioner(tntNetwork O11, tntNetwork O22, tntNetwork O12, tntNetwork O21, double tol){
        
    /* First declare required variables */
    tntNode beta, gamma, O, Osuper, invOsuper;
    tntNode OL11, OL22, OL12, OL21;
    tntNode row1, row2;
    int chiL1, chiR1, chiL2, chiR2;
    
    /* Set internal dimensions and check that internal dimensions of both fields 1 and 2 are equal */
    beta = tntNodeFindFirst(O11);
    O = tntNodeFindConn(beta,"S");
    gamma = tntNodeFindConn(O,"R");
    chiL1 = tntNodeGetLegDim(beta,"R");
    chiR1 = tntNodeGetLegDim(gamma, "L" );
    
    beta = tntNodeFindFirst(O22);
    O = tntNodeFindConn(beta,"S");
    gamma = tntNodeFindConn(O,"R");
    chiL2 = tntNodeGetLegDim(beta,"R");
    chiR2 = tntNodeGetLegDim(gamma, "L" );
    
    if(chiL1 != chiL2 || chiR1 != chiR2){
        printf("\nError: tntNodeNsIncomp2DConjGrad_BuildLeftEnvPreconditioner currently only works when internal dimensions of both fields are equal. Exiting.");
        exit(1);
    }
    
    /*
     |-----------------------------------------------------------------|
     | Steps:                                                          |
     | 1. Contract-out the right environment.                          |
     | 2. Direct sum O11 and O22 and anti-direct sum O12 and O21 to    |
     | produce Osuper := (OL11,OL12; OL22, OL21).                      |
     | 3. Invert Osuper with pseduo-inverse.                           |
     | 4. Free and exit.                                               |
     |-----------------------------------------------------------------|
    */
    
    /* Step 1: */
    OL11 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O11);
    OL12 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O12);
    OL21 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O21);
    OL22 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O22);
    
    /* Step 2: */
    tntNodeAddLeg(OL11,"ab"); tntNodeAddLeg(OL22,"ab");
    tntNodeAddLeg(OL12,"ab"); tntNodeAddLeg(OL21,"ab");
    
    row1 = tntNodeDirectSum(OL11,OL12,"a");
    row2 = tntNodeDirectSum(OL21,OL22,"a");
    Osuper = tntNodeDirectSum(row1, row2, "b");
    
    /* Step 3: Invert. Note: Osuper is destroyed by this operation. */
    invOsuper = tntNodePseudoInverse(Osuper,"TDb","RUa","RUaTDb=TDbRUa",tol); //Use TNT_DEFAULT_CG_RELTOL instead of TNT_DEFAULT_PSEUDOINVERSE_RELTOL.
    
    /* Step 4: Done. Now free everything, return the inverse and exit */
    tntNodeFree(&row1);
    tntNodeFree(&row2);
    tntNodeFree(&OL11);
    tntNodeFree(&OL22);
    tntNodeFree(&OL12);
    tntNodeFree(&OL21);
    
    return invOsuper;
}

/*! \ingroup ns
 * This function, along with the two others in this .c file, builds a preconditioner for use in the conjugate gradient descent optimisation algorithm of tntNodeNsIncomp2DConjGrad(...). Note: I have personally not  observed any significant speedup when using this preconditioner. Indeed, in my experience, it just slows the code down and therefore I reccommend against activating this preconditioner in tntNodeNsIncomp2DConjGrad(...).
 */

tntNode tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(tntNetwork Oorig){
        
    /* First declare required variables */
    tntNetwork Ocopy;
    tntNode V, beta, gamma, O;
    tntNode OL;
    
    /* Copy original variational networks */
    Ocopy = tntNetworkCopy(Oorig);

    /* Set the constituting nodes of the variational network */
    beta = tntNodeFindFirst(Ocopy);
    O = tntNodeFindConn(beta,"S");
    gamma = tntNodeFindConn(O,"R");
    
    /* Convert into node-group to make tensor algebra operations possible */
    tntNetworkToNodeGroup(&Ocopy, 0);
    
    /* Get rid of the singleton legs on the beta and gamma */
    tntNodeSqueeze(beta,"LMN");
    tntNodeSqueeze(gamma,"RST");
    /* Remove V from the network */
    V = tntNodeFindConn(O,"U");
    tntNodeRemove(V); tntNodeFree(&V);
    
    /* Contract-out the environment legs of the gammas */
    /* NOTE: Assumes legs "L" and "N" have same dimension. */
    tntNodeJoin(gamma,"L",gamma,"N");
    tntNodeContractSelf(&gamma);
    
    /* Contract the variational network into a node */
    OL = tntNodeListContract("RTDU",beta, O, gamma);
    
    /* Now finished. Nothing to free, so just return OL and exit. */
    return OL;
}

/*! \ingroup ns
 * This function, along with the two others in this .c file, builds a preconditioner for use in the conjugate gradient descent optimisation algorithm of tntNodeNsIncomp2DConjGrad(...). Note: I have personally not  observed any significant speedup when using this preconditioner. Indeed, in my experience, it just slows the code down and therefore I reccommend against activating this preconditioner in tntNodeNsIncomp2DConjGrad(...).
 */

tntNode tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(tntNetwork Oorig){
        
    /* First declare required variables */
    tntNetwork Ocopy;
    tntNode V, beta, gamma, O;
    tntNode OR;
    
    /* Copy original variational networks */
    Ocopy = tntNetworkCopy(Oorig);

    /* Set the constituting nodes of the variational network */
    beta = tntNodeFindFirst(Ocopy);
    O = tntNodeFindConn(beta,"S");
    gamma = tntNodeFindConn(O,"R");
    
    /* Convert into node-group to make tensor algebra operations possible */
    tntNetworkToNodeGroup(&Ocopy, 0);
    
    /* Get rid of the singleton legs on the beta and gamma */
    tntNodeSqueeze(beta,"LMN");
    tntNodeSqueeze(gamma,"RST");
    /* Remove V from the network */
    V = tntNodeFindConn(O,"U");
    tntNodeRemove(V); tntNodeFree(&V);
    
    /* Contract-out the environment legs of the betas */
    /* NOTE: Assumes legs "L" and "N" have same dimension. */
    tntNodeJoin(beta,"R",beta,"T");
    tntNodeContractSelf(&beta);
    
    /* Contract the variational network into a node */
    OR = tntNodeListContract("DULN",beta, O, gamma);
    
    /* Now finished. Nothing to free, so just return OR and exit. */
    return OR;
}
