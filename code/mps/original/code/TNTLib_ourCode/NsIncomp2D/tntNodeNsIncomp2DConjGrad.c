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
 * Solves linear problems of form M*X = B for X. Assumes M is a symmetric, positive semi-definite operator node. Assumes X and B are vector nodes. Assumes row legs of X and B are LRDa. Assumes row legs of M are MSDa, and col legs LRUb. "B" is provided as A1 and A2, to be direct-summed into Asuper = B, X as V1 and V2 from the networks, to be direct-summed into Vsuper = X, and M in the form of the identity plus the grand operator (O11, O12; O21, O22). The 4 sub-operators of this "grand operator" are taken from the inputed subnetworks of O11, O22, O12, O21.
 * Will modify all the inputed tntNode and tntNetwork objects.
 */

tntNode tntNodeNsIncomp2DConjGrad(tntNode A1, tntNode A2, tntNetwork O11,tntNetwork O22,tntNetwork O12,tntNetwork O21, int maxIter, double relTol, char preconditioner[]){
    
    /* Declare relevant variables */
    tntNode Vsuper, Vsuperc, Asuper;
    tntNode dummy, beta, V1, V2, V1_up, V2_up;
    tntComplexArray dummyArray;
    tntNode inversePreonditioner, inversePreonditionerc;
    
    struct tnsr_CG commObj; /* structure containg the results and necessary data for the reverse communication interface */
    TNT_ERR_RET_DEFS
    
    /* First extract V1 and V2 */
    beta = tntNodeFindFirst(O11);
    V1 = tntNodeCopy(tntNodeFindConn(beta,"R"));
    
    beta = tntNodeFindFirst(O22);
    V2 = tntNodeCopy(tntNodeFindConn(beta,"R"));
    
    /* Construct Vsuper and Asuper – the variational vector and the RHS vector. Use tntNodeReshape to keep their legs in a consistent order. */
    /* First Asuper: direct sum them to create the super-vector */
    tntNodeAddLeg(A1,"a"); tntNodeAddLeg(A2,"a");
    Asuper = tntNodeDirectSum(A1, A2, "a");
    
    /* Ditto for Vsuper*/
    tntNodeAddLeg(V1,"a"); tntNodeAddLeg(V2,"a");
    Vsuper = tntNodeDirectSum(V1, V2, "a");
    tntNodeFree(&V1); tntNodeFree(&V2);
    
    /* ----------------------------------------------- */
    /* Initialise CG method along with preconditioner  */
    /* ----------------------------------------------- */
    tntNodeReshapeAsMatrix(Vsuper,"DLR","a"); tntNodeReshapeAsMatrix(Asuper,"DLR","a"); //Make sure leg orders are consistent
        
    if( !strcmp(preconditioner,"Left Block-Inverse") ) {
        
        commObj.precond = 1; // Tell the routine to run the preconditioned variant of the algorithm.
        inversePreonditioner = tntNodeNsIncomp2DConjGrad_BuildLeftEnvPreconditioner(O11, O22, O12, O21, relTol);
    }
    else if( !strcmp(preconditioner,"None") ) {
        
        commObj.precond = 0; // Tell the routine to run the non-preconditioned variant of the algorithm.
    }
    else{
        printf("\nError: Preconditioning setting not recognised. Either use ""Left Block-Inverse"" or ""None"". ""Block Jacobi"" has been discontinued. Exiting. \n");
        exit(1);
    }
    
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
    dummyArray = tntComplexArrayAlloc(2,1);
    
    /* As long as action = 1 or 3, keep iterating. */
    while(commObj.action == 1 || commObj.action == 3){
        
        if(commObj.action == 1){
            
            /* Decompose Vsuper into V1 and V2 */
            Vsuperc = tntNodeCopy(Vsuper);
            
            /* First V1 */
            dummyArray.vals[1].re = 0;
            dummyArray.vals[0].re = 1;
            dummy = tntNodeCreate(&dummyArray,"b",2);
            tntNodeJoin(Vsuperc,"a",dummy,"b");
            V1 = tntNodeListContract("DLR",Vsuperc,dummy);
            
            /* Then get V2 */
            dummyArray.vals[0].re = 0;
            dummyArray.vals[1].re = 1;
            dummy = tntNodeCreate(&dummyArray,"b",2);
            tntNodeJoin(Vsuper,"a",dummy,"b");
            V2 = tntNodeListContract("DLR",Vsuper,dummy);
            
            /* Carry out contraction with super-operator */
            /* First contract O11*V1 and O22*V2  */
            V1_up = tntMpsVarMinMpo1sContract(V1, O11);
            V2_up = tntMpsVarMinMpo1sContract(V2, O22);
            
            /* Then apply the identity */
            tntNodeAdd(V1_up, V1);
            tntNodeAdd(V2_up, V2);
            
            /* And finally add O12*V2 and O21*V1 */
            dummy = tntMpsVarMinMpo1sContract(V2, O12);
            tntNodeAdd(V1_up, dummy); tntNodeFree(&dummy);
            dummy = tntMpsVarMinMpo1sContract(V1, O21);
            tntNodeAdd(V2_up, dummy); tntNodeFree(&dummy);
            
            /* V1 and V2 no longer needed */
            tntNodeFree(&V1); tntNodeFree(&V2);
            
            /* Build up Vsuper again */
            tntNodeAddLeg(V1_up,"a"); tntNodeAddLeg(V2_up,"a");
            Vsuper = tntNodeDirectSum(V1_up, V2_up, "a");
            
            /* V1_up and V2_up no longer needed */
            tntNodeFree(&V1_up); tntNodeFree(&V2_up);
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
            else{
                printf("\nError: Preconditioning setting not recognised. Either use ""Block Jacobi"" (NOT recommended), ""Left Block-Inverse"" or ""None"". Exiting. \n");
                exit(1);
            }
            
        }
        
        /* Done. Now feed the result back into the iteration routine. */
        tntNodeReshapeAsMatrix(Vsuper,"DLR","a"); //Set the indices of Vsuper to be in the same order as those of Asuper
        ret = _tnt_tnode_CG_iterate(Vsuper, Asuper, &commObj);
        TNT_PUB_ERR_CHK
        
    }
    
    /* This final Vsuper is the converged solution to the problem. Free everything and return Vsuper. */
    tntNodeFree(&Asuper);
    tntComplexArrayFree(&dummyArray);
    
    if( !strcmp(preconditioner,"Left Block-Inverse") ){
        tntNodeFree(&inversePreonditioner);
    }
    
    return Vsuper;
    
}
