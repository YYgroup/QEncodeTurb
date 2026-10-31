/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Need to configure tolerances for the local update algorithm */
#define TNT_DEFAULT_PSEUDOINVERSE_RELTOL 1e-8 /*! Default relative tolerance for setting singular values to zero during pseudoinverse. */
#define TNT_DEFAULT_CG_MAXITER 100 /*! Default number of permitted iterations for conjugate gradients algorithm. If preconditioner is used, 5 times more iterations will be allowed to justify its cost.*/
#define TNT_DEFAULT_CG_RELTOL 1e-8 /*! Default convergance criteria for conjugate gradients algorithm */

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 *
 * This function optimises the pair of tensors V1 and V2 (corresponding to flow along the x and y directions, respectively) for the overall optimisation scheme of tntMpsNsIncomp2D1sSweep(...). This is done by solving a local linear problem.
 *
 * \return No return value.
 */

void tntNodeNsIncomp2D1sLocalUpdate(tntNode *V1, tntNode *V2, tntNode U1, tntNode U2,tntNode UM1, tntNode UM2, tntNode *Bsum, tntNode D1xD1, tntNode D1xD2, tntNode D2xD2, tntNode R1xu1_left, tntNode R1xu1_right, tntNode *R1xBxu1_leftSum, tntNode *R1xBxu1_rightSum, tntNode R1d1xd1R1_left, tntNode R1d1xd1R1_right, tntNode R1d1xd2R2_left, tntNode R1d1xd2R2_right, tntNode R2xu2_left, tntNode R2xu2_right, tntNode *R2xBxu2_leftSum, tntNode *R2xBxu2_rightSum, tntNode R2d2xd2R2_left, tntNode R2d2xd2R2_right, unsigned numBurgTerms, char method[], char preconditioner[]){

    /* First declare relevant variables */
    tntNode V1_updated, V2_updated;
    tntNode A1, A2;
    tntNode A, O, beta, gamma;
    tntNode dummy;
    tntNode Vsuper, Vsuperc, Asuper;
    tntComplexArray dummyArray;
    
    /*
     |---------------------------------------------------------------|
     | Start by building up the vectors A1 and A2,                   |
     | and multiplying lambda*h^2 into D1xD1,D1xD2,D2xD2.            |
     | Will afterwards be used to solve the following coupled system |
     | of matrix equations, for V1 and V2:                           |
     | (I+lambda*O1)*V1 = A1 + lambda*P1*V2,                         |
     | (I+lambda*O2)*V2 = A2 + lambda*P2*V1.                         |
     | Note: D2xD1 and all related nodes are taken from the D1xD2    |
     | nodes. Simply exploit the fact that [d1,d2] = 0.              |
     |---------------------------------------------------------------|
    */
    
    /* Build up A1 and A2 */
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
    
    /* Now start contracting the Burger terms and adding them to A1 */
    for(unsigned i = 0; i<numBurgTerms; i++){
        
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
    
    /* Now start contracting the Burger terms and adding them to A2 */
    for(unsigned i = 0; i<numBurgTerms; i++){
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
    
    /* Now choose whether to use "inverse" (extremely inefficent since it scales as chi^6, but good for validation)
     or "Conjugate Gradient Descent" (quite efficient because it scales as chi^3). */
    if(!strcmp(method, "Inverse")){
        /* Solve by (pseudo) inverting an operator O: */
        
        /* First declare required variables */
        tntNode O1, P1, O2, P2, I;
        tntNode Osuper, Psuper;
        
        /* Build up O1 and O2 */
        /* First O1 */
        beta = tntNodeCopy(R1d1xd1R1_left);
        gamma = tntNodeCopy(R1d1xd1R1_right);
        /* Get rid of the singleton legs on beta and gamma */
        tntNodeSqueeze(beta,"LMN");
        tntNodeSqueeze(gamma,"RST");
        
        O = tntNodeCopy(D1xD1);
        tntNodeJoin(beta, "S", O, "L");
        tntNodeJoin(gamma, "M", O, "R");
        O1 = tntNodeListContract("LMRSDU", beta, gamma, O);
        
        /* Then O2 */
        beta = tntNodeCopy(R2d2xd2R2_left);
        gamma = tntNodeCopy(R2d2xd2R2_right);
        /* Get rid of the singleton legs on beta and gamma */
        tntNodeSqueeze(beta,"LMN");
        tntNodeSqueeze(gamma,"RST");
        
        O = tntNodeCopy(D2xD2);
        tntNodeJoin(beta, "S", O, "L");
        tntNodeJoin(gamma, "M", O, "R");
        O2 = tntNodeListContract("LMRSDU", beta, gamma, O);
        
        /* Build up P1 and P2 */
        /* First P1 */
        beta = tntNodeCopy(R1d1xd2R2_left);
        gamma = tntNodeCopy(R1d1xd2R2_right);
        /* Get rid of the singleton legs on beta and gamma */
        tntNodeSqueeze(beta,"LMN");
        tntNodeSqueeze(gamma,"RST");
        
        O = tntNodeCopy(D1xD2);
        tntNodeJoin(beta, "S", O, "L");
        tntNodeJoin(gamma, "M", O, "R");
        P1 = tntNodeListContract("LMRSDU", beta, gamma, O);
        
        /* And finally P2 */
        P2 = tntNodeCopy(P1);
        tntNodeMapLegs(P2,"LRUMSD=MSDLRU");
        
        /*
         |------------------------------------------------------------------|
         | Time to solve by use of pseudo-inverse              |
         | (I+lambda*O1)*V1 = A1+lambda*P1*V2,                              |
         | (I+lambda*O2)*V2 = A2+lambda*P2*V1.                              |
         | This is equivalent to solving                                    |
         | ( (I,0; 0,I) + lambda*(O1,-P1; -P2,O2) )*(V1; V2) = (A1; A2)     |
         |                                                                  |
         | =>                                                               |
         |                                                                  |
         | Steps:                                                           |
         | I0. Let I+O1 -> O1 and I+O2 -> O2                                |
         | I1. Direct sum A1 and A2 to produce Asuper = (A1; A2).           |
         | I2. Direct sum O1 and O2 and anti-direct sum P1 and P2 to produce|
         | Osuper := (O1,P1; P2, O2).                                       |
         | I3. Invert Osuper and multiply with (A1,A2) to give (V1,V2)      |
         | I4. Free and exit.                                               |
         |------------------------------------------------------------------|
        */
        
        /* Step I0: Create I and add to O1 and O2 */
        I = tntMpsCreateEye2(O1,"MSD","LRU");
        tntNodeAdd(O1,I);
        tntNodeAdd(O2,I);
        
        /* Step I1: */
        tntNodeAddLeg(A1,"a"); tntNodeAddLeg(A2,"a");
        Asuper = tntNodeDirectSum(A1, A2, "a");
        
        /* Step I2: */
        tntNodeAddLeg(O1,"ab"); tntNodeAddLeg(O2,"ab");
        tntNodeAddLeg(P1,"ab"); tntNodeAddLeg(P2,"ab");
        Osuper = tntNodeDirectSum(O1, O2, "ab");
        Psuper = tntNodeAntiDirectSum(P1,P2,"ab");
        tntNodeAdd(Osuper,Psuper);
        
        /* Step I3: Invert and compute Vsuper */
        tntNode OsuperInv = tntNodePseudoInverse(Osuper,"MSDa","LRUb","LRUbMSDa=MSDaLRUb",TNT_DEFAULT_PSEUDOINVERSE_RELTOL);
        tntNodeJoin(OsuperInv,"LRUb",Asuper,"LRDa");
        Vsuper = tntNodeListContract("LRDa", OsuperInv, Asuper); //Nik Gourianov note 01072020: replaced LRUa with LRDa.
        
        /* Step I4: Now free everything and exit */
        tntNodeFree(&I);
        tntNodeFree(&Psuper);
        tntNodeFree(&O1);
        tntNodeFree(&O2);
        tntNodeFree(&P1);
        tntNodeFree(&P2);
        
    }
    else if(!strcmp(method, "CG")){
        /* Solve iteratively with conjugate gradient descent. This is the reccommended approach,
         as it takes advantage of the symmetric & positive-semidefinite nautre of O. */
                
        /* First declare relevant variables */
        tntNetwork O1V1, O2V2, P1V2, P2V1;
       
        /* Produce networks which represent the operators O1, O2, P1, P2. */
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

        /* Now P1 */
        beta = tntNodeCopy(R1d1xd2R2_left);
        gamma = tntNodeCopy(R1d1xd2R2_right);
        /* Produce the network to be contracted: P1 with V2. */
        P1V2 = tntMpsVarMinMpo1sBuild(beta,gamma,*V2,D1xD2);
        /* And free beta & gamma before proceeding: */
        tntNodeFree(&beta); tntNodeFree(&gamma);
        
        /* And finally P2 */
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
        
        /* Feed everything into the conjugate gradient function */
        Vsuper = tntNodeNsIncomp2DConjGrad(A1, A2, O1V1, O2V2, P1V2, P2V1, TNT_DEFAULT_CG_MAXITER, TNT_DEFAULT_CG_RELTOL, preconditioner);
        
        /* Free the networks */
        tntNetworkFree(&O1V1);
        tntNetworkFree(&O2V2);
        tntNetworkFree(&P1V2);
        tntNetworkFree(&P2V1);
    }
    
    /* Extract V1 and V2 from Vsuper */
    /* First get V1 */
    Vsuperc = tntNodeCopy(Vsuper);
    dummyArray = tntComplexArrayAlloc(2,1);
    dummyArray.vals[0].re = 1;
    dummy = tntNodeCreate(&dummyArray,"b",2);
    
    tntNodeJoin(Vsuperc,"a",dummy,"b");
    V1_updated = tntNodeListContract("DLR",Vsuperc,dummy);
    
    /* Then get V2 */
    dummyArray.vals[0].re = 0;
    dummyArray.vals[1].re = 1;
    dummy = tntNodeCreate(&dummyArray,"b",2);
    
    tntNodeJoin(Vsuper,"a",dummy,"b");
    V2_updated = tntNodeListContract("DLR",Vsuper,dummy);
    
    /* Now replace V1 and V2 with V1_updated and V2_updated in their respective networks */
    tntNodeReplace(*V1,V1_updated);
    tntNodeReplace(*V2,V2_updated);
    tntNodeFree(V1);
    tntNodeFree(V2);
    *V1 = V1_updated;
    *V2 = V2_updated;
    
    /* Free everything which remains yet unfreed */
    tntNodeFree(&A1);
    tntNodeFree(&A2);
    tntComplexArrayFree(&dummyArray);
}
