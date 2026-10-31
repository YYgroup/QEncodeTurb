/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for NS and PDE functions, as well as headers for some internal functions of the TNT library */
#include "tntPde.h"
#include "tntIncompNs.h"
#include "../lib/headers/tnt_int.h"
#include "../lib/headers/dec_pde_int.h"

/*! \ingroup ns
 *
 * Solves linear problems of form M*X = B for X. Assumes M is a symmetric, positive semi-definite operator node. Assumes X and B are vector nodes. Assumes row legs of X and B are LRDa. Assumes row legs of M are MSDa, and col legs LRUb. "B" is provided as A1, A2 and A3, to be direct-summed into Asuper = B, X as V1, V2 and V3 from the networks to be direct-summed into Vsuper = X, and M in the form of the identity plus the grand operator (O1, P1V2, P1V3; P2V1, O2, P2V3; P3V1, P3V2, O3). The 9 sub-operators of this "grand operator" are taken from the inputed subnetworks of O1, O2, O3, P1V2, P1V3, P2V1, P2V3, P3V1, P3V2.
 * Will modify all the inputed tntNode and tntNetwork objects.
 */

tntNode tntNodeNsIncomp3DConjGrad(tntNode A1, tntNode A2,tntNode A3, tntNetwork O1, tntNetwork O2, tntNetwork O3, tntNetwork P1V2, tntNetwork P1V3,tntNetwork P2V1, tntNetwork P2V3,tntNetwork P3V1,tntNetwork P3V2, int maxIter, double relTol, char preconditioner[]){
    
    /* Declare relevant variables */
    tntNode Vsuper, Vsuperc, Asuper;
    tntNode dummy, beta, V1, V2, V3, V1_up, V2_up, V3_up;
    tntNode inversePreonditioner, inversePreonditionerc;
    tntComplexArray dummyArray;
    struct tnsr_CG commObj; /* structure containg the results and necessary data for the reverse communication interface */
    TNT_ERR_RET_DEFS
    
    /* First extract V1, V2 and V3 */
    beta = tntNodeFindFirst(O1);
    V1 = tntNodeCopy(tntNodeFindConn(beta,"R"));
    
    beta = tntNodeFindFirst(O2);
    V2 = tntNodeCopy(tntNodeFindConn(beta,"R"));
    
    beta = tntNodeFindFirst(O3);
    V3 = tntNodeCopy(tntNodeFindConn(beta,"R"));
    
    /* Construct Vsuper and Asuper – the variational vector and the RHS vector. Use tntNodeReshape to keep their legs in a consistent order. */
    /* First Asuper: direct sum them to create the super-vector */
    tntNodeAddLeg(A1,"a"); tntNodeAddLeg(A2,"a"); tntNodeAddLeg(A3,"a");
    dummy = tntNodeDirectSum(A1, A2, "a");
    Asuper = tntNodeDirectSum(dummy, A3, "a"); tntNodeFree(&dummy);
    
    /* Ditto for Vsuper*/
    tntNodeAddLeg(V1,"a"); tntNodeAddLeg(V2,"a"); tntNodeAddLeg(V3,"a");
    dummy = tntNodeDirectSum(V1, V2, "a");
    Vsuper = tntNodeDirectSum(dummy, V3, "a"); tntNodeFree(&dummy);
    tntNodeFree(&V1); tntNodeFree(&V2); tntNodeFree(&V3);
    
    /* ----------------------------------------------- */
    /* Initialise CG method along with preconditioner  */
    /* ----------------------------------------------- */
    
    if( !strcmp(preconditioner,"Left Block-Inverse") ) {
        
        commObj.precond = 1; // Tell the routine to run the preconditioned variant of the algorithm.
        inversePreonditioner = tntNodeNsIncomp3DConjGrad_BuildLeftEnvPreconditioner(O1, P1V2, P1V3, P2V1, O2 ,P2V3, P3V1, P3V2, O3, 1e-12);

    }
    else if( !strcmp(preconditioner,"Right Block-Inverse") ) {
        
        commObj.precond = 1; // Tell the routine to run the preconditioned variant of the algorithm.
        inversePreonditioner = tntNodeNsIncomp3DConjGrad_BuildRightEnvPreconditioner(O1, P1V2, P1V3, P2V1, O2 ,P2V3, P3V1, P3V2, O3, 1e-12);

    }
    else if( !strcmp(preconditioner,"None") ) {
        
        commObj.precond = 0; // Tell the routine to run the non-preconditioned variant of the algorithm.
    }
    else{
        printf("\nError: Preconditioning setting not recognised. Either use ""Left Block-Inverse"" or ""None"". ""Block Jacobi"" has been discontinued. Exiting. \n");
        exit(1);
    }
    
    /* Make sure leg orders are consistent */
    tntNodeReshapeAsMatrix(Vsuper,"DLR","a"); tntNodeReshapeAsMatrix(Asuper,"DLR","a");
    
    /* Now initialise the iterator */
    ret = _tnt_tnode_CG_init(maxIter, relTol, Vsuper, Asuper, &commObj);
    TNT_PUB_ERR_CHK
    
    /* And send the original vector into it */
    ret = _tnt_tnode_CG_iterate(Vsuper, Asuper, &commObj);
    TNT_PUB_ERR_CHK
    
    /* -------------------- */
    /* Now start iterating. */
    /* -------------------- */
    
    /* Allocate the array for picking out flow directions */
    dummyArray = tntComplexArrayAlloc(3,1);
    
    /* As long as action = 1 or 3, keep iterating. */
    while(commObj.action == 1 || commObj.action == 3){
        
        if(commObj.action == 1){
            
            /* Decompose Vsuper into V1 and V2 */
            Vsuperc = tntNodeCopy(Vsuper);
            /* First V1 */
            dummyArray.vals[2].re = 0;
            dummyArray.vals[0].re = 1;
            dummy = tntNodeCreate(&dummyArray,"b",3);
            tntNodeJoin(Vsuperc,"a",dummy,"b");
            V1 = tntNodeListContract("DLR",Vsuperc,dummy);
            
            /* Then get V2 */
            Vsuperc = tntNodeCopy(Vsuper);
            dummyArray.vals[0].re = 0;
            dummyArray.vals[1].re = 1;
            dummy = tntNodeCreate(&dummyArray,"b",3);
            tntNodeJoin(Vsuperc,"a",dummy,"b");
            V2 = tntNodeListContract("DLR",Vsuperc,dummy);
            
            /* And finally V3 */
            dummyArray.vals[1].re = 0;
            dummyArray.vals[2].re = 1;
            dummy = tntNodeCreate(&dummyArray,"b",3);
            tntNodeJoin(Vsuper,"a",dummy,"b");
            V3 = tntNodeListContract("DLR",Vsuper,dummy);
            
            /* Carry out contraction */
            /* First contract O1V1, O2V2 and O3V3 */
            V1_up = tntMpsVarMinMpo1sContract(V1, O1);
            V2_up = tntMpsVarMinMpo1sContract(V2, O2);
            V3_up = tntMpsVarMinMpo1sContract(V3, O3);
            
            /* Then apply the identity */
            tntNodeAdd(V1_up, V1);
            tntNodeAdd(V2_up, V2);
            tntNodeAdd(V3_up, V3);
            
            /* And finally add P1V2, P2V1 etc.*/
            /* First for V1 */
            dummy = tntMpsVarMinMpo1sContract(V2, P1V2);
            tntNodeAdd(V1_up, dummy); tntNodeFree(&dummy);
            dummy = tntMpsVarMinMpo1sContract(V3, P1V3);
            tntNodeAdd(V1_up, dummy); tntNodeFree(&dummy);
            
            /* Then for V2 */
            dummy = tntMpsVarMinMpo1sContract(V1, P2V1);
            tntNodeAdd(V2_up, dummy); tntNodeFree(&dummy);
            dummy = tntMpsVarMinMpo1sContract(V3, P2V3);
            tntNodeAdd(V2_up, dummy); tntNodeFree(&dummy);
            
            /* And finally for V3 */
            dummy = tntMpsVarMinMpo1sContract(V1, P3V1);
            tntNodeAdd(V3_up, dummy); tntNodeFree(&dummy);
            dummy = tntMpsVarMinMpo1sContract(V2, P3V2);
            tntNodeAdd(V3_up, dummy); tntNodeFree(&dummy);
            
            /* V1, V2 and V3 no longer needed */
            tntNodeFree(&V1); tntNodeFree(&V2); tntNodeFree(&V3);
            
            /* Build up Vsuper again */
            tntNodeAddLeg(V1_up,"a"); tntNodeAddLeg(V2_up,"a"); tntNodeAddLeg(V3_up,"a");
            dummy = tntNodeDirectSum(V1_up, V2_up, "a");
            Vsuper = tntNodeDirectSum(dummy, V3_up, "a"); tntNodeFree(&dummy);
            
            /* V1_up, V2_up and V3_up no longer needed */
            tntNodeFree(&V1_up); tntNodeFree(&V2_up); tntNodeFree(&V3_up);
            
        }
        /* In case preconditioning is turned on, make use of it. The code should only trigger if "(commObj.ipar)[10] = 1".  */
        else if(commObj.action == 3){
            
            if( !strcmp(preconditioner,"Left Block-Inverse") ) {
                inversePreonditionerc = tntNodeCopy(inversePreonditioner);
                Vsuperc = tntNodeCopy(Vsuper); tntNodeFree(&Vsuper);
                tntNodeJoin(inversePreonditionerc,"RUa", Vsuperc,"LDa");
                Vsuper = tntNodeListContract("DLaR", inversePreonditionerc, Vsuperc);
                
                //Note, Nik: Strange way of doing it, I know, but have to due to TNT's habit of shallow-copying of tnodes.
            }
            else if( !strcmp(preconditioner,"Right Block-Inverse") ) {
                inversePreonditionerc = tntNodeCopy(inversePreonditioner);
                Vsuperc = tntNodeCopy(Vsuper); tntNodeFree(&Vsuper);
                tntNodeJoin(inversePreonditionerc,"LUa", Vsuperc,"RDa");
                Vsuper = tntNodeListContract("DRaL", inversePreonditionerc, Vsuperc);
                
                //Note, Nik: Strange way of doing it, I know, but have to due to TNT's habit of shallow-copying of tnodes.
            }
            else{
                printf("\nError: Preconditioning setting not recognised. Either use ""Block Jacobi"" (NOT recommended), ""Left Block-Inverse"" or ""None"". Exiting. \n");
                exit(1);
            }
            
        }
        
        /* Done. Now feed the result back into the iteration routine. */
        tntNodeReshapeAsMatrix(Vsuper,"DLR","a");
        ret = _tnt_tnode_CG_iterate(Vsuper, Asuper, &commObj);
        TNT_PUB_ERR_CHK
    }
    
    /* This final Vsuper is the converged solution to the problem. Free everything and return Vsuper. */
    tntNodeFree(&Asuper);
    if( !strcmp(preconditioner,"Left Block-Inverse")  || !strcmp(preconditioner,"Right Block-Inverse") ){
        tntNodeFree(&inversePreonditioner);
    }

    tntComplexArrayFree(&dummyArray);
    
    return Vsuper;
    
}
