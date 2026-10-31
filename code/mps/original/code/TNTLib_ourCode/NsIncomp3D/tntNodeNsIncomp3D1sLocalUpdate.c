/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Need to configure tolerances for the local update algorithm */
#define TNT_DEFAULT_CG_MAXITER 150 /*! Maximum number of permitted iterations for conjugate gradients algorithm */
#define TNT_DEFAULT_CG_RELTOL 1e-8 /*! Default convergance criteria for conjugate gradients algorithm */

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 *
 * This function optimises the pair of tensors V1, V2 and V3 (corresponding to flow along the x, y and z directions, respectively) for the overall optimisation scheme of tntMpsNsIncomp3D1sSweep(...). This is done by solving a local linear problem.
 *
 * \return No return value.
 */

void tntNodeNsIncomp3D1sLocalUpdate(tntNode *V1, tntNode *V2, tntNode *V3, tntNode U1, tntNode U2, tntNode U3, tntNode UM1, tntNode UM2, tntNode UM3, tntNode *Bsum,
                                    tntNode D1xD1, tntNode D1xD2, tntNode D1xD3, tntNode D2xD2, tntNode D2xD3, tntNode D3xD3,
                                    tntNode R1xu1_left, tntNode R1xu1_right, tntNode *R1xBxu1_leftSum, tntNode *R1xBxu1_rightSum, tntNode R1d1xd1R1_left, tntNode R1d1xd1R1_right, tntNode R1d1xd2R2_left, tntNode R1d1xd2R2_right, tntNode R1d1xd3R3_left, tntNode R1d1xd3R3_right,
                                    tntNode R2xu2_left, tntNode R2xu2_right, tntNode *R2xBxu2_leftSum, tntNode *R2xBxu2_rightSum, tntNode R2d2xd2R2_left, tntNode R2d2xd2R2_right, tntNode R2d2xd3R3_left, tntNode R2d2xd3R3_right,
                                    tntNode R3xu3_left, tntNode R3xu3_right, tntNode *R3xBxu3_leftSum, tntNode *R3xBxu3_rightSum, tntNode R3d3xd3R3_left, tntNode R3d3xd3R3_right,
                                    unsigned numExTerms, char preconditioner[]){
    
    /* First declare relevant variables */
    tntNode V1_updated, V2_updated, V3_updated;
    tntNode A1, A2, A3;
    tntNode A, O, beta, gamma;
    tntNode dummy;
    tntNode Vsuper, Vsuperc, Asuper;
    tntComplexArray dummyArray;
    
/*
     |---------------------------------------------------------------------------|
     | Willsolve the following coupled system of equations                       |
     | of matrix equations, for V1, V2 and V3:                                   |
     | (I+lambda*O1)*V1 = A1 + lambda*(P1V2+P1V3),                               |
     | (I+lambda*O2)*V2 = A2 + lambda*(P2V1+P3V3),                               |
     | (I+lambda*O3)*V3 = A3 + lambda*(P3V1+P3V2).                               |
     | Assumes lambda is already pre-multiplied into D1xD1, D1xD2, etc.          |
     | Note: D2xD1 and all related nodes are taken from the D1xD2 nodes.         |
     | Ditto with the D3xD1 and D3xD2 nodes. Simply exploit the fact that        |
     | [d1,d2], [d1,d3] and [d2,d3] are all equal to 0.                          |
     |---------------------------------------------------------------------------|
*/
    
    /* Build up A1, A2 and A3 */
    /* Start with A1 */
    /* Initialise A1 with the contraction of R1xu1_left, R1xu1_right and U1 */
    beta = tntNodeCopy(R1xu1_left);
    gamma = tntNodeCopy(R1xu1_right);
    /* Get rid of the singleton legs on beta and gamma */
    tntNodeSqueeze(beta,"LM");
    tntNodeSqueeze(gamma,"RS");
    
    A = tntNodeCopy(U1);
    tntNodeJoin(beta, "R", A, "L");
    tntNodeJoin(gamma, "L", A, "R");
    A1 = tntNodeListContract("LRD", beta, gamma, A);
    
    /* Now start contracting the advection terms and adding them to A1 */
    for(unsigned i = 0; i<numExTerms; i++){
        
        beta = tntNodeCopy(R1xBxu1_leftSum[i]);
        gamma = tntNodeCopy(R1xBxu1_rightSum[i]);
        /* Get rid of the singleton legs on beta and gamma */
        tntNodeSqueeze(beta,"LMN");
        tntNodeSqueeze(gamma,"RST");
        
        A = tntNodeCopy(UM1);
        O = tntNodeCopy(Bsum[i]);
        tntNodeJoin(beta, "R", A, "L");
        tntNodeJoin(gamma, "L", A, "R");
        tntNodeJoin(beta, "S", O, "L");
        tntNodeJoin(gamma, "M", O, "R");
        tntNodeJoin(A, "D", O, "U");
        
        dummy = tntNodeListContract("LRD", beta, gamma, A, O);
        tntNodeAdd(A1, dummy);
        tntNodeFree(&dummy);        
    }
    
    /* Now build A2 */
    /* Initialise A2 with the contraction of R2xu2_left, R2xu2_right and U2 */
    beta = tntNodeCopy(R2xu2_left);
    gamma = tntNodeCopy(R2xu2_right);
    /* Get rid of the singleton legs on beta and gamma */
    tntNodeSqueeze(beta,"LM");
    tntNodeSqueeze(gamma,"RS");
    
    A = tntNodeCopy(U2);
    tntNodeJoin(beta, "R", A, "L");
    tntNodeJoin(gamma, "L", A, "R");
    A2 = tntNodeListContract("LRD", beta, gamma, A);
    
    /* Now start contracting the advection terms and adding them to A2 */
    for(unsigned i = 0; i<numExTerms; i++){
        
        beta = tntNodeCopy(R2xBxu2_leftSum[i]);
        gamma = tntNodeCopy(R2xBxu2_rightSum[i]);
        /* Get rid of the singleton legs on beta and gamma */
        tntNodeSqueeze(beta,"LMN");
        tntNodeSqueeze(gamma,"RST");
        
        A = tntNodeCopy(UM2);
        O = tntNodeCopy(Bsum[i]);
        tntNodeJoin(beta, "R", A, "L");
        tntNodeJoin(gamma, "L", A, "R");
        tntNodeJoin(beta, "S", O, "L");
        tntNodeJoin(gamma, "M", O, "R");
        tntNodeJoin(A, "D", O, "U");
        
        dummy = tntNodeListContract("LRD", beta, gamma, A, O);
        tntNodeAdd(A2, dummy);
        tntNodeFree(&dummy);
    }
    
    /* And finally do A3 */
    /* Initialise A3 with the contraction of R3xu3_left, R3xu3_right and U3 */
    beta = tntNodeCopy(R3xu3_left);
    gamma = tntNodeCopy(R3xu3_right);
    /* Get rid of the singleton legs on beta and gamma */
    tntNodeSqueeze(beta,"LM");
    tntNodeSqueeze(gamma,"RS");
    
    A = tntNodeCopy(U3);
    tntNodeJoin(beta, "R", A, "L");
    tntNodeJoin(gamma, "L", A, "R");
    A3 = tntNodeListContract("LRD", beta, gamma, A);
    
    /* Now start contracting the advection terms and adding them to A3 */
    for(unsigned i = 0; i<numExTerms; i++){
        
        beta = tntNodeCopy(R3xBxu3_leftSum[i]);
        gamma = tntNodeCopy(R3xBxu3_rightSum[i]);
        /* Get rid of the singleton legs on beta and gamma */
        tntNodeSqueeze(beta,"LMN");
        tntNodeSqueeze(gamma,"RST");
        
        A = tntNodeCopy(UM3);
        O = tntNodeCopy(Bsum[i]);
        tntNodeJoin(beta, "R", A, "L");
        tntNodeJoin(gamma, "L", A, "R");
        tntNodeJoin(beta, "S", O, "L");
        tntNodeJoin(gamma, "M", O, "R");
        tntNodeJoin(A, "D", O, "U");
        
        dummy = tntNodeListContract("LRD", beta, gamma, A, O);
        tntNodeAdd(A3, dummy);
        tntNodeFree(&dummy);
    }
    /* Solve iteratively with conjugate gradient descent. This is a good approach,
     as it takes advantage of the symmetric & positive-semidefinite nautre of O. */
    
    /* First declare relevant variables */
    tntNetwork O1V1, O2V2, O3V3, P1V2, P1V3, P2V1, P2V3, P3V1, P3V2;
    
    /* Produce networks which represent the operators */
    /* Start with O1 */
    beta = tntNodeCopy(R1d1xd1R1_left);
    gamma = tntNodeCopy(R1d1xd1R1_right);
    /* Produce the network to be contracted: O1 with V1. */
    O1V1 = tntMpsVarMinMpo1sBuild(beta,gamma,*V1,D1xD1);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now O2 */
    beta = tntNodeCopy(R2d2xd2R2_left);
    gamma = tntNodeCopy(R2d2xd2R2_right);
    /* Produce the network to be contracted: O2 with V2. */
    O2V2 = tntMpsVarMinMpo1sBuild(beta,gamma,*V2,D2xD2);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now O3 */
    beta = tntNodeCopy(R3d3xd3R3_left);
    gamma = tntNodeCopy(R3d3xd3R3_right);
    /* Produce the network to be contracted: O2 with V2. */
   // tntNodePrintAll(*V3);
    O3V3 = tntMpsVarMinMpo1sBuild(beta,gamma,*V3,D3xD3);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now P1V2 */
    beta = tntNodeCopy(R1d1xd2R2_left);
    gamma = tntNodeCopy(R1d1xd2R2_right);
    /* Produce the network to be contracted: P1 with V2. */
    P1V2 = tntMpsVarMinMpo1sBuild(beta,gamma,*V2,D1xD2);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now P1V3 */
    beta = tntNodeCopy(R1d1xd3R3_left);
    gamma = tntNodeCopy(R1d1xd3R3_right);
    /* Produce the network to be contracted: P1 with V3. */
    P1V3 = tntMpsVarMinMpo1sBuild(beta,gamma,*V3,D1xD3);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now P2V1 */
    beta = tntNodeCopy(R1d1xd2R2_left);
    gamma = tntNodeCopy(R1d1xd2R2_right);
    O = tntNodeCopy(D1xD2);
    tntNodeMapLegs(beta,"LMNRST=NMLTSR");
    tntNodeMapLegs(gamma,"RSTLMN=TSRNML");
    tntNodeMapLegs(O,"LRUD=LRDU");
    /* Produce the network to be contracted: P2 with V1. */
    P2V1 = tntMpsVarMinMpo1sBuild(beta,gamma,*V1,O);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now P2V3 */
    beta = tntNodeCopy(R2d2xd3R3_left);
    gamma = tntNodeCopy(R2d2xd3R3_right);
    O = tntNodeCopy(D2xD3);
    /* Produce the network to be contracted: P2 with V3. */
    P2V3 = tntMpsVarMinMpo1sBuild(beta,gamma,*V3,O);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Now P3V1 */
    beta = tntNodeCopy(R1d1xd3R3_left);
    gamma = tntNodeCopy(R1d1xd3R3_right);
    O = tntNodeCopy(D1xD3);
    tntNodeMapLegs(beta,"LMNRST=NMLTSR");
    tntNodeMapLegs(gamma,"RSTLMN=TSRNML");
    tntNodeMapLegs(O,"LRUD=LRDU");
    /* Produce the network to be contracted: P2 with V1. */
    P3V1 = tntMpsVarMinMpo1sBuild(beta,gamma,*V1,O);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* And finally, P3V2 */
    beta = tntNodeCopy(R2d2xd3R3_left);
    gamma = tntNodeCopy(R2d2xd3R3_right);
    O = tntNodeCopy(D2xD3);
    tntNodeMapLegs(beta,"LMNRST=NMLTSR");
    tntNodeMapLegs(gamma,"RSTLMN=TSRNML");
    tntNodeMapLegs(O,"LRUD=LRDU");
    /* Produce the network to be contracted: P2 with V3. */
    P3V2 = tntMpsVarMinMpo1sBuild(beta,gamma,*V2,O);
    /* And free beta & gamma before proceeding: */
    tntNodeFree(&beta); tntNodeFree(&gamma);
    
    /* Feed everything into the conjugate gradient function */
    Vsuper = tntNodeNsIncomp3DConjGrad(A1, A2, A3, O1V1, O2V2, O3V3, P1V2, P1V3, P2V1, P2V3, P3V1, P3V2, TNT_DEFAULT_CG_MAXITER, TNT_DEFAULT_CG_RELTOL, preconditioner); //Vsuper has legs DLRa.
        
    /* Free the networks */
    tntNetworkFree(&O1V1);
    tntNetworkFree(&O2V2);
    tntNetworkFree(&O3V3);
    tntNetworkFree(&P1V2);
    tntNetworkFree(&P1V3);
    tntNetworkFree(&P2V1);
    tntNetworkFree(&P2V3);
    tntNetworkFree(&P3V1);
    tntNetworkFree(&P3V2);
    
    /* Extract V1 and V2 from Vsuper */
    /* First get V1 */
    Vsuperc = tntNodeCopy(Vsuper);
    dummyArray = tntComplexArrayAlloc(3,1);
    dummyArray.vals[0].re = 1;
    dummy = tntNodeCreate(&dummyArray,"b",3);

    tntNodeJoin(Vsuperc,"a",dummy,"b");
    V1_updated = tntNodeListContract("DLR",Vsuperc,dummy);

    /* Then get V2 */
    Vsuperc = tntNodeCopy(Vsuper);
    dummyArray.vals[0].re = 0;
    dummyArray.vals[1].re = 1;
    dummy = tntNodeCreate(&dummyArray,"b",3);
    
    tntNodeJoin(Vsuperc,"a",dummy,"b");
    V2_updated = tntNodeListContract("DLR",Vsuperc,dummy);
    
    /* Then get V3 */
    dummyArray.vals[1].re = 0;
    dummyArray.vals[2].re = 1;
    dummy = tntNodeCreate(&dummyArray,"b",3);
    
    tntNodeJoin(Vsuper,"a",dummy,"b");
    V3_updated = tntNodeListContract("DLR",Vsuper,dummy);
    
    /* Now replace V1, V2 and V3 with V1_updated, V2_updated and V3_updated in their respective networks */
    tntNodeReplace(*V1,V1_updated);
    tntNodeReplace(*V2,V2_updated);
    tntNodeReplace(*V3,V3_updated);
    tntNodeFree(V1);
    tntNodeFree(V2);
    tntNodeFree(V3);
    *V1 = V1_updated;
    *V2 = V2_updated;
    *V3 = V3_updated;

    /* Free everything which remains yet unfreed */
    tntNodeFree(&A1);
    tntNodeFree(&A2);
    tntNodeFree(&A3);
    tntComplexArrayFree(&dummyArray);
    
}
