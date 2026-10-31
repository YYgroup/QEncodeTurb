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
 * This function, along with the other one in this .c file, builds a preconditioner for use in the conjugate gradient descent optimisation algorithm of tntNodeNsIncomp3DConjGrad(...). Note: I have personally not  observed any significant speedup when using this preconditioner. Indeed, in my experience, it just slows the code down and therefore I reccommend against activating this preconditioner in tntNodeNsIncomp3DConjGrad(...).
 */

tntNode tntNodeNsIncomp3DConjGrad_BuildLeftEnvPreconditioner(tntNetwork O11, tntNetwork O12, tntNetwork O13, tntNetwork O21, tntNetwork O22, tntNetwork O23, tntNetwork O31, tntNetwork O32, tntNetwork O33, double tol){
        
    /* First declare required variables */
    tntNode dummy, beta, gamma, O, Osuper, invOsuper;
    tntNode zero, row1, row2, row3;
    tntNode OL11, OL12, OL13, OL21, OL22, OL23, OL31, OL32, OL33;
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
        printf("\nError: tntNodeNsIncomp3DConjGrad_BuildLeftEnvPreconditioner currently only works when internal dimensions of both fields are equal. Exiting.");
        exit(1);
    }

    /*
     |-----------------------------------------------------------------|
     | Steps:                                                          |
     | 1. Contract-out the right environment.                          |
     | 2. Direct sum O11, O22 and O33, and anti-direct sum O12, O21,   |
     | O13, O31, O23, O32 to produce                                   |
     | Osuper :=(OL11, OL12, OL13; OL21, OL22, OL23; OL31, OL32, OL33).|
     | 3. Invert Osuper with pseduo-inverse.                           |
     | 4. Free and exit.                                               |
     |-----------------------------------------------------------------|
    */

    /* Step 1: */
    OL11 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O11);
    OL12 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O12);
    OL13 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O13);
    OL21 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O21);
    OL22 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O22);
    OL23 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O23);
    OL31 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O31);
    OL32 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O32);
    OL33 = tntNodeNsIncompConjGrad_BuildLeftEnvPreconditioner_ExtractLeftEnv(O33);

    /* Step 2: Create 3 sub block-nodes*/
    /* Step 2a: Upper-left block */
    tntNodeAddLeg(OL11,"ab"); tntNodeAddLeg(OL12,"ab"); tntNodeAddLeg(OL13,"ab");
    tntNodeAddLeg(OL21,"ab"); tntNodeAddLeg(OL22,"ab"); tntNodeAddLeg(OL23,"ab");
    tntNodeAddLeg(OL31,"ab"); tntNodeAddLeg(OL32,"ab"); tntNodeAddLeg(OL33,"ab");
    
    dummy = tntNodeDirectSum(OL11,OL12,"a");
    row1 = tntNodeDirectSum(dummy,OL13,"a"); tntNodeFree(&dummy);
    
    dummy = tntNodeDirectSum(OL21,OL22,"a");
    row2 = tntNodeDirectSum(dummy,OL23,"a"); tntNodeFree(&dummy);
    
    dummy = tntNodeDirectSum(OL31,OL32,"a");
    row3 = tntNodeDirectSum(dummy,OL33,"a"); tntNodeFree(&dummy);
    
    dummy = tntNodeDirectSum(row1,row2,  "b");
    Osuper = tntNodeDirectSum(dummy,row3,"b"); tntNodeFree(&dummy);
    
    /* Step 3: Invert. Note: Osuper is destroyed by this operation. */
    invOsuper = tntNodePseudoInverse(Osuper,"TDb","RUa","RUaTDb=TDbRUa",tol); //Use TNT_DEFAULT_CG_RELTOL instead of TNT_DEFAULT_PSEUDOINVERSE_RELTOL.

    /* Step 4: Done. Now free everything, return the inverse and exit */
    tntNodeFree(&row1);
    tntNodeFree(&row2);
    tntNodeFree(&row3);
    tntNodeFree(&OL11);
    tntNodeFree(&OL12);
    tntNodeFree(&OL13);
    tntNodeFree(&OL21);
    tntNodeFree(&OL22);
    tntNodeFree(&OL23);
    tntNodeFree(&OL31);
    tntNodeFree(&OL32);
    tntNodeFree(&OL33);

    return invOsuper;
}

/*! \ingroup ns
 * This function, along with the other one in this .c file, builds a preconditioner for use in the conjugate gradient descent optimisation algorithm of tntNodeNsIncomp3DConjGrad(...). Note: I have personally not  observed any significant speedup when using this preconditioner. Indeed, in my experience, it just slows the code down and therefore I reccommend against activating this preconditioner in tntNodeNsIncomp3DConjGrad(...).
 */

tntNode tntNodeNsIncomp3DConjGrad_BuildRightEnvPreconditioner(tntNetwork O11, tntNetwork O12, tntNetwork O13, tntNetwork O21, tntNetwork O22, tntNetwork O23, tntNetwork O31, tntNetwork O32, tntNetwork O33, double tol){
        
    /* First declare required variables */
    tntNode dummy, beta, gamma, O, Osuper, invOsuper;
    tntNode zero, row1, row2, row3;
    tntNode OR11, OR12, OR13, OR21, OR22, OR23, OR31, OR32, OR33;
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
        printf("\nError: tntNodeNsIncomp3DConjGrad_BuildRightEnvPreconditioner currently only works when internal dimensions of both fields are equal. Exiting.");
        exit(1);
    }

    /*
     |-----------------------------------------------------------------|
     | Steps:                                                          |
     | 1. Contract-out the right environment.                          |
     | 2. Direct sum O11, O22 and O33, and anti-direct sum O12, O21,   |
     | O13, O31, O23, O32 to produce                                   |
     | Osuper :=(OR11, OR12, OR13; OR21, OR22, OR23; OR31, OR32, OR33).|
     | 3. Invert Osuper with pseduo-inverse.                           |
     | 4. Free and exit.                                               |
     |-----------------------------------------------------------------|
    */

    /* Step 1: */
    OR11 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O11);
    OR12 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O12);
    OR13 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O13);
    OR21 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O21);
    OR22 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O22);
    OR23 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O23);
    OR31 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O31);
    OR32 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O32);
    OR33 = tntNodeNsIncompConjGrad_BuildRightEnvPreconditioner_ExtractRightEnv(O33);

    /* Step 2: Create 3 sub block-nodes*/
    /* Step 2a: Upper-left block */
    tntNodeAddLeg(OR11,"ab"); tntNodeAddLeg(OR12,"ab"); tntNodeAddLeg(OR13,"ab");
    tntNodeAddLeg(OR21,"ab"); tntNodeAddLeg(OR22,"ab"); tntNodeAddLeg(OR23,"ab");
    tntNodeAddLeg(OR31,"ab"); tntNodeAddLeg(OR32,"ab"); tntNodeAddLeg(OR33,"ab");
    
    dummy = tntNodeDirectSum(OR11,OR12,"a");
    row1 = tntNodeDirectSum(dummy,OR13,"a"); tntNodeFree(&dummy);
    
    dummy = tntNodeDirectSum(OR21,OR22,"a");
    row2 = tntNodeDirectSum(dummy,OR23,"a"); tntNodeFree(&dummy);
    
    dummy = tntNodeDirectSum(OR31,OR32,"a");
    row3 = tntNodeDirectSum(dummy,OR33,"a"); tntNodeFree(&dummy);
    
    dummy = tntNodeDirectSum(row1,row2,  "b");
    Osuper = tntNodeDirectSum(dummy,row3,"b"); tntNodeFree(&dummy);
    
    /* Step 3: Invert. Note: Osuper is destroyed by this operation. */
    invOsuper = tntNodePseudoInverse(Osuper,"NDb","LUa","LUaNDb=NDbLUa",tol); //Use TNT_DEFAULT_CG_RELTOL instead of TNT_DEFAULT_PSEUDOINVERSE_RELTOL.

    /* Step 4: Done. Now free everything, return the inverse and exit */
    tntNodeFree(&row1);
    tntNodeFree(&row2);
    tntNodeFree(&row3);
    tntNodeFree(&OR11);
    tntNodeFree(&OR12);
    tntNodeFree(&OR13);
    tntNodeFree(&OR21);
    tntNodeFree(&OR22);
    tntNodeFree(&OR23);
    tntNodeFree(&OR31);
    tntNodeFree(&OR32);
    tntNodeFree(&OR33);

    return invOsuper;
}
