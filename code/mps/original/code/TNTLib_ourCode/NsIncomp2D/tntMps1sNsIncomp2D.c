/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for the PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 * Performs a single, optimising sweep of the MPS (explicit time-stepping). *
 * \return No return value.
 */
void tntMpsNsIncomp2D1sSweep(tntNetwork v1, /*!< Assumes canonical centre at site "start_site" */
                             tntNetwork v2, /*!< Assumes canonical centre at site "start_site" */
                             tntNetwork u1, tntNetwork u2, tntNetwork uMidpoint1,tntNetwork uMidpoint2, tntNetwork* BurgMpoSum,
                             tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d2xd2,
                             tntNodeArray *R1xu1_left, tntNodeArray *R1xu1_right,
                             tntNodeArray **R1xBxumid1_leftSum, tntNodeArray **R1xBxumid1_rightSum,
                             tntNodeArray *R2xu2_left, tntNodeArray *R2xu2_right,
                             tntNodeArray **R2xBxumid2_leftSum, tntNodeArray **R2xBxumid2_rightSum,
                             tntNodeArray *R1d1xd1R1_left, tntNodeArray *R1d1xd1R1_right,
                             tntNodeArray *R1d1xd2R2_left, tntNodeArray *R1d1xd2R2_right,
                             tntNodeArray *R2d2xd2R2_left, tntNodeArray *R2d2xd2R2_right,
                             tntLegLabel dir,
                             unsigned numBurgTerms,
                             int start_site,
                             int end_site,
                             char method[],
                             char preconditioner[]){
 
    
    /* Declare relevant variables */
    int L, j; /* Number of sites and site index. Here j is chosen to always refer to the left-most node in the pair of nodes being updated. */
    tntNode V1, V2; //Nodes to be time-evolved
    tntNode U1, U2, UM1, UM2;// Used to update L-1 node
    tntNode Q, R; //Variables used to move orthogonal centre
    
    /* These variables are needed to send in the sum of advection-diffusion terms into the local tntNode update function tntIncomp2D1sLocalTE */
    tntNode *Bsum;
    tntNode *R1xBxumid1_leftSumValsAtj, *R1xBxumid1_rightSumValsAtj;
    tntNode *R2xBxumid2_leftSumValsAtj, *R2xBxumid2_rightSumValsAtj;
    tntNode D1xD1, D1xD2, D2xD2;
    
    L = tntMpsLength(v1);
    
    /* Allocate memory to the node arrays */
    Bsum = malloc(numBurgTerms*sizeof(tntNode));
    R1xBxumid1_leftSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R1xBxumid1_rightSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R2xBxumid2_leftSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R2xBxumid2_rightSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));

    if ('R' == dir[0]) {
        /* Set the current tensors to be the first of their respective networks, except the advection-diffusion terms */
        V1 = tntNodeFindFirst(v1);
        V2 = tntNodeFindFirst(v2);
        U1 = tntNodeFindFirst(u1);
        U2 = tntNodeFindFirst(u2);
        UM1 = tntNodeFindFirst(uMidpoint1);
        UM2 = tntNodeFindFirst(uMidpoint2);
        D1xD1 = tntNodeFindFirst(d1xd1);
        D1xD2 = tntNodeFindFirst(d1xd2);
        D2xD2 = tntNodeFindFirst(d2xd2);
        
        /* Move forward until current tensors are start tensors, except the advection-diffusion terms */
        for (j = 0; j < start_site; j++) {
            V1 = tntNodeFindConn(V1,"R");
            V2 = tntNodeFindConn(V2,"R");
            U1 = tntNodeFindConn(U1,"R");
            U2 = tntNodeFindConn(U2,"R");
            UM1 = tntNodeFindConn(UM1,"R");
            UM2 = tntNodeFindConn(UM2,"R");
            
            D1xD1 = tntNodeFindConn(D1xD1,"R");
            D1xD2 = tntNodeFindConn(D1xD2,"R");
            D2xD2 = tntNodeFindConn(D2xD2,"R");
        }
        /* Advection-diffusion terms need to be treated specially as the advection-diffusion tensor network is represented as a direct sum of terms (and this direct sum is not taken on the network level for reasons of speed) */
        /* First set them to be the first of their respective networks */
        for(unsigned i = 0; i<numBurgTerms; i++){
            Bsum[i] = tntNodeFindFirst(BurgMpoSum[i]);
        }
        
        /* Then move forward */
        for (j = 0; j < start_site; j++) {
            for(unsigned i = 0; i<numBurgTerms; i++){
                Bsum[i] = tntNodeFindConn(Bsum[i],"R");
            }
        }
        
        for (j = start_site; j <= end_site-1; j++) {
            
            /* Gather the precontracted advection-diffusion networks into arrays before time-evolving */
            for(unsigned i = 0; i<numBurgTerms; i++){
                R1xBxumid1_leftSumValsAtj[i] = ((*R1xBxumid1_leftSum)[i]).vals[j];
                R1xBxumid1_rightSumValsAtj[i] = ((*R1xBxumid1_rightSum)[i]).vals[j];
                R2xBxumid2_leftSumValsAtj[i] = ((*R2xBxumid2_leftSum)[i]).vals[j];
                R2xBxumid2_rightSumValsAtj[i] = ((*R2xBxumid2_rightSum)[i]).vals[j];
            }
 
            
            /* Then perform a 1-site local update */
            tntNodeNsIncomp2D1sLocalUpdate(&V1, &V2,
                                       U1, U2, UM1, UM2, Bsum, D1xD1, D1xD2, D2xD2,
                                       R1xu1_left->vals[j], R1xu1_right->vals[j], R1xBxumid1_leftSumValsAtj, R1xBxumid1_rightSumValsAtj,
                                       R1d1xd1R1_left->vals[j], R1d1xd1R1_right->vals[j], R1d1xd2R2_left->vals[j], R1d1xd2R2_right->vals[j],
                                       R2xu2_left->vals[j], R2xu2_right->vals[j], R2xBxumid2_leftSumValsAtj, R2xBxumid2_rightSumValsAtj,
                                       R2d2xd2R2_left->vals[j], R2d2xd2R2_right->vals[j],
                                       numBurgTerms,method, preconditioner);

            
            /* Shift canonical centre from j to j+1 of field 1 */
            Q = tntNodeQR(V1,"LD", "R", "L");
            R = tntNodeFindConn(Q,"R");
            V1 = tntNodeContract(R, tntNodeFindConn(R,"R")); //Canonical centre now at site j+1
            V1 = tntNodeFindConn(V1,"L"); //Go back to site j
            
            /* Shift canonical centre from j to j+1 of field 2 */
            Q = tntNodeQR(V2,"LD", "R", "L");
            R = tntNodeFindConn(Q,"R");
            V2 = tntNodeContract(R, tntNodeFindConn(R,"R")); //Canonical centre now at site j+1
            V2 = tntNodeFindConn(V2,"L"); //Go back to site j
            
            /* Update node arrays */
            /* Start with updating the MPO ones */
            tntMpsPdeMpoUpdate(V1, D1xD1, V1, "R", j, R1d1xd1R1_left);
            tntMpsPdeMpoUpdate(V2, D1xD2, V1, "R", j, R1d1xd2R2_left);
            tntMpsPdeMpoUpdate(V2, D2xD2, V2, "R", j, R2d2xd2R2_left);
            
            /* Then update the MPS ones between v and u fields (between updated and original fields) */
            tntMpsPdeMpsUpdate(U1, V1, "R", j, R1xu1_left);
            tntMpsPdeMpsUpdate(U2, V2, "R", j, R2xu2_left);
            
            /* Finally, update the MPS ones between v and the advection-diffusion terms */
            for(unsigned i = 0; i<numBurgTerms; i++){
                tntMpsPdeMpoUpdate(UM1, Bsum[i], V1, "R", j, &((*R1xBxumid1_leftSum)[i]) );
                tntMpsPdeMpoUpdate(UM2, Bsum[i], V2, "R", j, &((*R2xBxumid2_leftSum)[i]) );
            }
            
            /* Go to next set of nodes for all the MPSes*/
            V1 = tntNodeFindConn(V1,"R");
            V2 = tntNodeFindConn(V2,"R");
            U1 = tntNodeFindConn(U1,"R");
            U2 = tntNodeFindConn(U2,"R");
            UM1 = tntNodeFindConn(UM1,"R");
            UM2 = tntNodeFindConn(UM2,"R");
            D1xD1 = tntNodeFindConn(D1xD1,"R");
            D1xD2 = tntNodeFindConn(D1xD2,"R");
            D2xD2 = tntNodeFindConn(D2xD2,"R");
            for(unsigned i = 0; i<numBurgTerms; i++){
                Bsum[i] = tntNodeFindConn(Bsum[i],"R");
            }
            
        }
    }
    else{
        /* Set the current tensors to be the last of their respective networks, except the advection-diffusion terms */
        V1 = tntNodeFindLast(v1);
        V2 = tntNodeFindLast(v2);
        U1 = tntNodeFindLast(u1);
        U2 = tntNodeFindLast(u2);
        UM1 = tntNodeFindLast(uMidpoint1);
        UM2 = tntNodeFindLast(uMidpoint2);
        D1xD1 = tntNodeFindLast(d1xd1);
        D1xD2 = tntNodeFindLast(d1xd2);
        D2xD2 = tntNodeFindLast(d2xd2);
        
        /* Move backwards until current tensors are start tensors, except the advection-diffusion terms */
        for (j = L-1; j > start_site; j--) {
            V1 = tntNodeFindConn(V1,"L");
            V2 = tntNodeFindConn(V2,"L");
            U1 = tntNodeFindConn(U1,"L");
            U2 = tntNodeFindConn(U2,"L");
            UM1 = tntNodeFindConn(UM1,"L");
            UM2 = tntNodeFindConn(UM2,"L");
            D1xD1 = tntNodeFindConn(D1xD1,"L");
            D1xD2 = tntNodeFindConn(D1xD2,"L");
            D2xD2 = tntNodeFindConn(D2xD2,"L");
        }
        /* Advection-diffusion terms need to be treated specially as the advection-diffusion tensor network is represented as a direct sum of terms (and this direct sum is not taken on the network level for reasons of speed) */
        for (unsigned i = 0; i<numBurgTerms; i++){
            /* First set them to be the last of their respective networks */
            Bsum[i] = tntNodeFindLast(BurgMpoSum[i]);
        }
        
        /* Then move backwards */
        for (j = L-1; j > start_site; j--){
            for (unsigned i = 0; i<numBurgTerms; i++){
                Bsum[i] = tntNodeFindConn(Bsum[i],"L");
            }
        }
        
        for (j = start_site; j >= end_site+1; j--) {
            
            /* Gather the precontracted advection-diffusion networks into arrays before time-evolving */
            for(unsigned i = 0; i<numBurgTerms; i++){
                R1xBxumid1_leftSumValsAtj[i] = ((*R1xBxumid1_leftSum)[i]).vals[j];
                R1xBxumid1_rightSumValsAtj[i] = ((*R1xBxumid1_rightSum)[i]).vals[j];
                R2xBxumid2_leftSumValsAtj[i] = ((*R2xBxumid2_leftSum)[i]).vals[j];
                R2xBxumid2_rightSumValsAtj[i] = ((*R2xBxumid2_rightSum)[i]).vals[j];
                
            }
            
            tntNodeNsIncomp2D1sLocalUpdate(&V1, &V2,
                                       U1, U2, UM1, UM2, Bsum, D1xD1, D1xD2, D2xD2,
                                       R1xu1_left->vals[j], R1xu1_right->vals[j], R1xBxumid1_leftSumValsAtj, R1xBxumid1_rightSumValsAtj,
                                       R1d1xd1R1_left->vals[j], R1d1xd1R1_right->vals[j], R1d1xd2R2_left->vals[j], R1d1xd2R2_right->vals[j],
                                       R2xu2_left->vals[j], R2xu2_right->vals[j], R2xBxumid2_leftSumValsAtj, R2xBxumid2_rightSumValsAtj,
                                       R2d2xd2R2_left->vals[j], R2d2xd2R2_right->vals[j],
                                       numBurgTerms, method, preconditioner);

            
            /* Shift canonical centre from j to j-1 of field 1 */
            Q = tntNodeQR(V1,"RD", "L", "R");
            R = tntNodeFindConn(Q,"L");
            V1 = tntNodeContract(R, tntNodeFindConn(R,"L")); //Canonical centre now at site j-1
            V1 = tntNodeFindConn(V1,"R"); //Go back to site j
            
            /* Shift canonical centre from j to j-1 of field 2 */
            Q = tntNodeQR(V2,"RD", "L", "R");
            R = tntNodeFindConn(Q,"L");
            V2 = tntNodeContract(R, tntNodeFindConn(R,"L")); //Canonical centre now at site j-1
            V2 = tntNodeFindConn(V2,"R"); //Go back to site j
            
            /* Update node arrays */
            /* Start with updating the MPO ones */
            tntMpsPdeMpoUpdate(V1, D1xD1, V1, "L", j, R1d1xd1R1_right);
            tntMpsPdeMpoUpdate(V2, D1xD2, V1, "L", j, R1d1xd2R2_right);
            tntMpsPdeMpoUpdate(V2, D2xD2, V2, "L", j, R2d2xd2R2_right);
            
            /* Then update the MPS ones between v and u fields (between updated and original fields) */
            tntMpsPdeMpsUpdate(U1,V1, "L", j, R1xu1_right);
            tntMpsPdeMpsUpdate(U2,V2, "L", j, R2xu2_right);
            
            /* Finally, upadte the MPS ones between v and the advection-diffusion terms */
            for(unsigned i = 0; i<numBurgTerms; i++){
                tntMpsPdeMpoUpdate(UM1, Bsum[i],V1, "L", j, &((*R1xBxumid1_rightSum)[i]));
                tntMpsPdeMpoUpdate(UM2, Bsum[i],V2, "L", j, &((*R2xBxumid2_rightSum)[i]));
            }
            
            /* Go to next set of nodes for all the MPSes*/
            V1 = tntNodeFindConn(V1,"L");
            V2 = tntNodeFindConn(V2,"L");
            U1 = tntNodeFindConn(U1,"L");
            U2 = tntNodeFindConn(U2,"L");
            UM1 = tntNodeFindConn(UM1,"L");
            UM2 = tntNodeFindConn(UM2,"L");
            D1xD1 = tntNodeFindConn(D1xD1,"L");
            D1xD2 = tntNodeFindConn(D1xD2,"L");
            D2xD2 = tntNodeFindConn(D2xD2,"L");
            for(unsigned i = 0; i<numBurgTerms; i++){
                Bsum[i] = tntNodeFindConn(Bsum[i],"L");
            }
        }
    }
    
    /* Must time-evolve final site as well. */
    /* Gather the precontracted advection-diffusion networks into arrays before time-evolving */
    for(unsigned i = 0; i<numBurgTerms; i++){
        R1xBxumid1_leftSumValsAtj[i] = ((*R1xBxumid1_leftSum)[i]).vals[end_site];
        R1xBxumid1_rightSumValsAtj[i] = ((*R1xBxumid1_rightSum)[i]).vals[end_site];
        R2xBxumid2_leftSumValsAtj[i] = ((*R2xBxumid2_leftSum)[i]).vals[end_site];
        R2xBxumid2_rightSumValsAtj[i] = ((*R2xBxumid2_rightSum)[i]).vals[end_site];
    }
    
    tntNodeNsIncomp2D1sLocalUpdate(&V1, &V2,
                               U1, U2,UM1, UM2, Bsum, D1xD1, D1xD2, D2xD2,
                               R1xu1_left->vals[end_site], R1xu1_right->vals[end_site], R1xBxumid1_leftSumValsAtj, R1xBxumid1_rightSumValsAtj,
                               R1d1xd1R1_left->vals[end_site], R1d1xd1R1_right->vals[end_site], R1d1xd2R2_left->vals[end_site], R1d1xd2R2_right->vals[end_site],
                               R2xu2_left->vals[end_site], R2xu2_right->vals[end_site], R2xBxumid2_leftSumValsAtj, R2xBxumid2_rightSumValsAtj,
                               R2d2xd2R2_left->vals[end_site], R2d2xd2R2_right->vals[end_site],
                               numBurgTerms,method, preconditioner);

    
    /* Free the arrays of nodes */
    free(Bsum);
    free(R1xBxumid1_leftSumValsAtj);
    free(R1xBxumid1_rightSumValsAtj);
    free(R2xBxumid2_leftSumValsAtj);
    free(R2xBxumid2_rightSumValsAtj);
}

/*! \ingroup ns
 * Performs a double-sweep optimising the MPS (explicit time-stepping).*
 * \return No return value.
 */

void tntMpsNsIncomp2D1sFullSweep(tntNetwork v1, tntNetwork v2,
                               tntNetwork u1, tntNetwork u2, tntNetwork uMidpoint1, tntNetwork uMidpoint2,
                               tntNetwork *BurgMpoSum,
                               tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d2xd2,
                               tntNodeArray *R1xu1_left, tntNodeArray *R1xu1_right,
                               tntNodeArray **R1xBxumid1_leftSum, tntNodeArray **R1xBxumid1_rightSum,
                               tntNodeArray *R2xu2_left, tntNodeArray *R2xu2_right,
                               tntNodeArray **R2xBxumid2_leftSum, tntNodeArray **R2xBxumid2_rightSum,
                               tntNodeArray *R1d1xd1R1_left, tntNodeArray *R1d1xd1R1_right,
                               tntNodeArray *R1d1xd2R2_left, tntNodeArray *R1d1xd2R2_right,
                               tntNodeArray *R2d2xd2R2_left, tntNodeArray *R2d2xd2R2_right,
                               unsigned numBurgTerms, /* Number of terms in the advection-diffusion term 1 & 2 sums */
                               char method[],
                               char preconditioner[]){
    
    unsigned L = tntMpsLength(v1);
    
    if(L == 1){
        tntErrorPrint("Cannot integrate the incompressible NS equation forward if it just has one site! Exiting."); /* NO_COVERAGE */
        exit(-1);
    }
    else{
        /* Initialise relevant nodes */
        tntNode Q, R, V1, V2; /* Needed to shift canonical centre from L-1 to L-2 */
        tntNode D1xD1, D1xD2, D2xD2;

        /* Sweep from node 0 to node L-1, i.e. from the very first to the very last */
        tntMpsNsIncomp2D1sSweep(v1, v2,
                                u1, u2, uMidpoint1, uMidpoint2, BurgMpoSum,
                                d1xd1, d1xd2, d2xd2,
                                R1xu1_left, R1xu1_right, R1xBxumid1_leftSum, R1xBxumid1_rightSum,
                                R2xu2_left, R2xu2_right, R2xBxumid2_leftSum, R2xBxumid2_rightSum,
                                R1d1xd1R1_left, R1d1xd1R1_right, R1d1xd2R2_left, R1d1xd2R2_right,
                                R2d2xd2R2_left, R2d2xd2R2_right,
                                "R", numBurgTerms, 0, L-1, method, preconditioner);

        /* Shift canonical centre from L-1 to L-2 of field 1 */
        V1 = tntNodeFindLast(v1);
        Q = tntNodeQR(V1,"RD", "L", "R");
        R = tntNodeFindConn(Q,"L");
        V1 = tntNodeContract(R, tntNodeFindConn(R,"L")); //Canonical centre now at site L-2
        
        /* Shift canonical centre from L-1 to L-2 of field 2 */
        V2 = tntNodeFindLast(v2);
        Q = tntNodeQR(V2,"RD", "L", "R");
        R = tntNodeFindConn(Q,"L");
        V2 = tntNodeContract(R, tntNodeFindConn(R,"L")); //Canonical centre now at site L-2

        /* Update node arrays */
        /* Go back to L-1 */
        V1 = tntNodeFindConn(V1,"R");
        V2 = tntNodeFindConn(V2,"R");
        
        /* Then do it for the sum of the advection-diffusion terms */
        for(unsigned i = 0; i<numBurgTerms; i++){
            tntMpsPdeMpoUpdate(tntNodeFindLast(uMidpoint1),tntNodeFindLast(BurgMpoSum[i]),V1, "L", L-1, &((*R1xBxumid1_rightSum)[i]) );
            tntMpsPdeMpoUpdate(tntNodeFindLast(uMidpoint2),tntNodeFindLast(BurgMpoSum[i]),V2, "L", L-1, &((*R2xBxumid2_rightSum)[i]) );
        }

        /* Then for the differencing networks: */
        D1xD1 = tntNodeFindLast(d1xd1);
        D1xD2 = tntNodeFindLast(d1xd2);
        D2xD2 = tntNodeFindLast(d2xd2);
        tntMpsPdeMpoUpdate(V1, D1xD1, V1, "L", L-1, R1d1xd1R1_right);
        tntMpsPdeMpoUpdate(V2, D1xD2, V1, "L", L-1, R1d1xd2R2_right);
        tntMpsPdeMpoUpdate(V2, D2xD2, V2, "L", L-1, R2d2xd2R2_right);

        /* And finally for the initial fields: */
        tntMpsPdeMpsUpdate(tntNodeFindLast(u1),V1, "L", L-1, R1xu1_right);
        tntMpsPdeMpsUpdate(tntNodeFindLast(u2),V2, "L", L-1, R2xu2_right);

        /* Now Sweep from node L-2 to node 0, i.e. from the second-to-last to the first, and finish */
        tntMpsNsIncomp2D1sSweep(v1, v2,
                                u1, u2, uMidpoint1, uMidpoint2, BurgMpoSum,
                                d1xd1, d1xd2, d2xd2,
                                R1xu1_left, R1xu1_right, R1xBxumid1_leftSum, R1xBxumid1_rightSum,
                                R2xu2_left, R2xu2_right, R2xBxumid2_leftSum, R2xBxumid2_rightSum,
                                R1d1xd1R1_left, R1d1xd1R1_right, R1d1xd2R2_left, R1d1xd2R2_right,
                                R2d2xd2R2_left, R2d2xd2R2_right,
                                "L", numBurgTerms, L-2, 0, method, preconditioner);
    }
    
}

/*! \ingroup ns
 * Performs full Runge-Kutta 2 time-step  of the MPS (explicit time-stepping).
 * \return No return value.
 */

void tntMpsNsIncomp2D1sTimeStep(tntNetwork v1, tntNetwork v2,
                            tntNetwork u1, tntNetwork u2, tntNetwork uMidpoint1, tntNetwork uMidpoint2,
                            tntNetwork *BurgMpoSum,
                            tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d2xd2,
                            tntNodeArray *R1d1xd1R1_left, tntNodeArray *R1d1xd1R1_right,
                            tntNodeArray *R1d1xd2R2_left, tntNodeArray *R1d1xd2R2_right,
                            tntNodeArray *R2d2xd2R2_left, tntNodeArray *R2d2xd2R2_right,
                            unsigned numBurgTerms, /* Number of terms in the Burger MPOs sum array */
                            double epsilon /* convergence number for sweeping */){

    /* Declare relevant variables */
    int numRuns;
    double Ek0, Ek1 = 0; //Kinetic energy of field
    unsigned L = tntMpsLength(v1);
    tntNodeArray R1xu1_left, R1xu1_right;
    tntNodeArray R2xu2_left, R2xu2_right;
    tntNodeArray *R1xBxumid1_leftSum, *R1xBxumid1_rightSum;
    tntNodeArray *R2xBxumid2_leftSum, *R2xBxumid2_rightSum;
    tntNetwork burgMpo;
    
    /* Allocate memory to the arrays of node arrays */
    R1xBxumid1_leftSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R1xBxumid1_rightSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R2xBxumid2_leftSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R2xBxumid2_rightSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    
    /* Produce the pre-contracted networks */
    tntMpsPdeMpsInit(u1,v1,&R1xu1_left,&R1xu1_right);
    tntMpsPdeMpsInit(u2,v2,&R2xu2_left,&R2xu2_right);
    
    for(unsigned i = 0; i<numBurgTerms; i++){
        burgMpo = BurgMpoSum[i];
        tntMpsPdeMpoInit(uMidpoint1,burgMpo,v1,&(R1xBxumid1_leftSum[i]),&(R1xBxumid1_rightSum[i]));
        tntMpsPdeMpoInit(uMidpoint2,burgMpo,v2,&(R2xBxumid2_leftSum[i]),&(R2xBxumid2_rightSum[i]));
    }
    
    /* Sweep until convergence */
    numRuns = 0;
    Ek0 = 0;
    Ek1 = 2*epsilon;
    while( fabs( ((Ek1-Ek0)/Ek0)) > epsilon || numRuns<2){

        tntMpsNsIncomp2D1sFullSweep(v1, v2,
                                    u1, u2, uMidpoint1, uMidpoint2,
                                    BurgMpoSum,
                                    d1xd1, d1xd2, d2xd2,
                                    &R1xu1_left, &R1xu1_right,
                                    &R1xBxumid1_leftSum, &R1xBxumid1_rightSum,
                                    &R2xu2_left, &R2xu2_right,
                                    &R2xBxumid2_leftSum, &R2xBxumid2_rightSum,
                                    R1d1xd1R1_left, R1d1xd1R1_right,
                                    R1d1xd2R2_left, R1d1xd2R2_right,
                                    R2d2xd2R2_left, R2d2xd2R2_right,
                                    numBurgTerms, "CG", "None");//"Left Block-Inverse");//  "None" , "Left Block-Inverse".
        numRuns++;
        
        /* Update Kinetic energy */
        Ek0=Ek1;
        Ek1 = tntMpsSelfProduct(v1,-1)+tntMpsSelfProduct(v2,-1);
        printf("\nDouble-Sweep completed. Ek1=%.6e, (Ek1-Ek0)/Ek0=%.6e", Ek1, (Ek1-Ek0)/Ek0 );
    }
    
    printf("\n=====Convergence reached at epsilon=%.2e after %i runs=====",epsilon, numRuns);
    
    /* Free the node arrays */
    tntNodeArrayFree(&R1xu1_left); tntNodeArrayFree(&R1xu1_right);
    tntNodeArrayFree(&R2xu2_left); tntNodeArrayFree(&R2xu2_right);
    
    for(unsigned i = 0; i<numBurgTerms; i++){
        tntNodeArrayFree( &(R1xBxumid1_leftSum[i]));
        tntNodeArrayFree( &(R1xBxumid1_rightSum[i]));
        tntNodeArrayFree( &(R2xBxumid2_leftSum[i]));
        tntNodeArrayFree( &(R2xBxumid2_rightSum[i]));
    }

    /* Free the array of node arrays */
    free(R1xBxumid1_leftSum); free(R1xBxumid1_rightSum);
    free(R2xBxumid2_leftSum); free(R2xBxumid2_rightSum);
    
}


/*! \ingroup ns
 * Performs multiple Runge-Kutta 2 time-steps of the MPS (explicit time-stepping).
 * \return No return value.
 */

void tntMpsBoxTurbulence2DTimeEvolutionRK2(tntNetwork *v1, tntNetwork *v2,/* Inputed MPS flow field  */
                                           double dt, /* Time step */
                                           double T, /* Final time */
                                           double T_sampling, /* Sampling period (i.e. 1/T_sampling = sampling frequency) */
                                           char savepath[TNT_STRLEN],  /* Output path */
                                           char BC1[], char BC2[], /* Boundary conditions of the box */
                                           double Re, /* Reynolds number */
                                           double lambda, /* Penalty coefficient */
                                           double epsilon, /* Convergence sensitivity */
                                           int chi, int chi_mpo, /* Bond dimension for MPSes and MPOs, respectively */
                                           double *errp /* Dummy variable */){
    
    extern tntSystem systype;
    
    
    /* First of all, make sure the user has turned of symmetries. */
    if(TNT_SYMM_NONE != tntSymmTypeGet()){
        char warning_msg[TNT_STRLEN];
        strcpy(warning_msg,"tntMpsBoxTurbulence2DTimeEvolutionRK2 doesn't work with symmetries! Turning them OFF and continuing!");
        tntSymmTypeUnset();
        tntWarningPrint(warning_msg);
    }
    
    /* Declare relevant variables */
    char currentSavename[TNT_STRLEN];
    char fieldname[TNT_STRLEN];
    tntNetwork u1, u2, uMidpoint1, uMidpoint2;
    tntNetwork *BurgMpoSum;
    tntNetwork d1, d2;
    tntNetwork dd1, dd2;
    tntNetwork d1xd1, d1xd2, d2xd2;
    tntNodeArray R1d1xd1R1_left, R1d1xd1R1_right;
    tntNodeArray R1d1xd2R2_left, R1d1xd2R2_right;
    tntNodeArray R2d2xd2R2_left, R2d2xd2R2_right;
    tntNode A, B;
    double chi_Lleg, chi_Rleg;
    double t;
    double lengthScale = 1;
    unsigned numBurgTerms;
    unsigned L = tntMpsLength(*v1);
    int step, stepsBetweenSampling, L_dim, R_dim;
    double h  = pow((double) 1/2,(double) L);
    
    /* Sampling variables */
    char txtfilename[TNT_STRLEN];
    tntDoubleArray Ekin, Enst;
    
    /* Contract the original field down to the appropriate bond dimension and use it as the "starting guess" for the time-evolved field. Also, put it in appropriate orthogonal form. */
    tntMpsOrth(*v1,L-1);
    tntMpsOrth(*v2,L-1);
    *errp = *errp + tntMpsTruncate(*v1, 0, chi);
    *errp = *errp + tntMpsTruncate(*v2, 0, chi);
    
    /* Initialise differencing MPO & operators*/
    d1 = tntMpsPde2DCreateOrder1Diff1(BC1,L);
    d2 = tntMpsPde2DCreateOrder1Diff2(BC2,L);
    dd1 = tntMpsPde2DCreateOrder2Diff1(BC1,L);
    dd2 = tntMpsPde2DCreateOrder2Diff2(BC2,L);
    
    /* Create MPO networks as products between differencers and their adjoints */
    d1xd1 = tntMpoMpoProduct(d1, d1);
    d2xd2 = tntMpoMpoProduct(d2, d2);
    d1xd2 = tntMpoMpoProduct(d2, d1);

    /* Exploit sparseness to decrease bond dimension */
    *errp += tntMpoTruncate(d1xd1, -1);
    *errp += tntMpoTruncate(d2xd2, -1);
    *errp += tntMpoTruncate(d1xd2, -1);
    
    /* Scale by -1 to account for the fact that the differencing operator is anti-symmetric, and multiply by lambda and by dx^2 so that lagrange multiplier in cost function has same units as the error */
    tntNodeScaleReal(tntNodeFindFirst(d1xd1),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d2xd2),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d1xd2),-lambda);
    
    /* Initialise Kinetic energy and Enstrophy arrays */
    Ekin = tntDoubleArrayAlloc(ceil(T/dt)+1);
    Enst = tntDoubleArrayAlloc(ceil(T/dt)+1);
    
    /* Prepare system for time-evolution */
    t=0; //Start at time = 0.
    step = 0;
    stepsBetweenSampling = floor(T_sampling/dt);
    if(stepsBetweenSampling == 0){
        stepsBetweenSampling = 1;
    }
    
    /* Initialise initial field */
    u1 = tntNetworkCopy(*v1);
    u2 = tntNetworkCopy(*v2);
    
    /* Set initial fields to be random to support prolongation */
    tntNetworkFree(v1);
    tntNetworkFree(v2);
    *v1 = tntNetworkCreate();
    *v2 = tntNetworkCreate();
    R_dim = 1;
    for(unsigned i = 0;i<L; i++){
        /* Set dimensions*/
        L_dim = R_dim;
        R_dim = MIN(MIN(pow(4,i+1), pow(4,L-i-1)), chi);

        /* Insert nodes */
        A = tntNodeCreate(NULL,"DLR",4,L_dim, R_dim);
        B = tntNodeCreate(NULL,"DLR",4,L_dim, R_dim);
        tntNodeInsertAtEnd(A,"L","R",*v1);
        tntNodeInsertAtEnd(B,"L","R",*v2);
    }
    tntMpsOrth(*v1,0);
    tntMpsOrth(*v2,0);

    /* Initialise pre-contracted networks */
    tntMpsPdeMpoInit(*v1,d1xd1,*v1,&R1d1xd1R1_left,&R1d1xd1R1_right);
    tntMpsPdeMpoInit(*v2,d1xd2,*v1,&R1d1xd2R2_left,&R1d1xd2R2_right);
    tntMpsPdeMpoInit(*v2,d2xd2,*v2,&R2d2xd2R2_left,&R2d2xd2R2_right);
    
    /* Save initial field */
    printf("\nSaving initial fields at t=%.0f. Ek0 = %f.",t,tntMpsSelfProduct(u1,0)+ tntMpsSelfProduct(u2,0)); fflush(stdout);
    currentSavename[0] = '\0';
    strcat(currentSavename,savepath);
    sprintf(fieldname, "/NsIncomp2DFieldOutL%iC%iRe%.0fmt%.0f.mat",L,chi,Re,round(t*1000));
    strcat(currentSavename, fieldname );
    tntNetworksSave(currentSavename, u1, u2);
    
    /* Calculate initial kinetic energy and enstrophy */
    tntMpsNsIncomp2DCalcEkinEnst(&(Ekin.vals[0]),&(Enst.vals[0]),u1,u2,d1,d2,chi,0);
    
    /* Run RK2, i.e. midpoint method, to time-evolve */
    time_t start, end;
    double time_used;
    while(t<T*0.99999 || T == 0){
        start = time(NULL);
        tntWalltimeSet(0.0);

        /* First calculate advection-diffusion term at initial field and store in ad1 and ad2 arrays */
        tntMpsNsIncomp2DBuildBurgerTerms(&BurgMpoSum, &numBurgTerms,
                                     u1, u2,
                                     d1, dd1,
                                     d2, dd2,
                                     Re, dt/2, chi,chi_mpo, errp);

        /* Now do a time-step dt/2 forward */
        tntMpsNsIncomp2D1sTimeStep(*v1, *v2,
                               u1, u2, u1, u2,
                               BurgMpoSum,
                               d1xd1, d1xd2, d2xd2,
                               &R1d1xd1R1_left, &R1d1xd1R1_right,
                               &R1d1xd2R2_left, &R1d1xd2R2_right,
                               &R2d2xd2R2_left, &R2d2xd2R2_right,
                               numBurgTerms,
                               epsilon);
        
        /* Extract midpoint field */
        uMidpoint1 = tntNetworkCopy(*v1);
        uMidpoint2 = tntNetworkCopy(*v2);
                
        /* Now delete old advection-diffusion terms and calculate new advection-diffusion term at the current field, i.e. at the midpoint */
        tntNetworkArrayFree(BurgMpoSum, numBurgTerms);
        tntMpsNsIncomp2DBuildBurgerTerms(&BurgMpoSum, &numBurgTerms,
                                         uMidpoint1, uMidpoint2,
                                         d1, dd1,
                                         d2, dd2,
                                         Re, dt, chi, chi_mpo, errp);

        /* Use the updated advection-diffusion terms to time-evolve a full dt forward */
        tntMpsNsIncomp2D1sTimeStep(*v1, *v2,
                               u1, u2, uMidpoint1, uMidpoint2,
                               BurgMpoSum,
                               d1xd1, d1xd2, d2xd2,
                               &R1d1xd1R1_left, &R1d1xd1R1_right,
                               &R1d1xd2R2_left, &R1d1xd2R2_right,
                               &R2d2xd2R2_left, &R2d2xd2R2_right,
                               numBurgTerms,
                               epsilon);
        
        
        /* Move one timestep forward */
        tntNetworkFree(&u1);
        tntNetworkFree(&u2);
        tntNetworkFree(&uMidpoint1);
        tntNetworkFree(&uMidpoint2);
        u1 = tntNetworkCopy(*v1);
        u2 = tntNetworkCopy(*v2);
        t=t+dt;
        step++;
        
        /* Advection-diffusion terms now out of date, so must be deleted in preperation for next step */
        tntNetworkArrayFree(BurgMpoSum, numBurgTerms);
        
        end = time(NULL);
        time_used = ((double) (end - start));// / CLOCKS_PER_SEC;
        printf("\nTotal time spent doing one RK2 timestep: %.0fs. Of that, %.0fs was spent in tensor_reshape.",time_used, tntWalltimeGet()); fflush(stdout);
        /* If dt = 0, the time-evolution is finished. Get out of the loop. */
        if(t==0 || T==0){
            break;
        }
        
        /* Calculate kinetic energy and enstrophy at time t. */
        tntMpsNsIncomp2DCalcEkinEnst(&(Ekin.vals[step]),&(Enst.vals[step]),u1,u2,d1,d2,chi,0);
        
        printf("\nRK2 step finished at t = %f, of T = %f. Ek1 = %.4f\n",t, T, tntMpsSelfProduct(u1,0)+ tntMpsSelfProduct(u2,0)); fflush(stdout);
        
        /* If appropriate, sample the current MPS and return */
        if(step%stepsBetweenSampling == 0){
            printf("Saving time-evolved fields at t=%f.\n",t);
            currentSavename[0] = '\0';
            strcat(currentSavename,savepath);
            sprintf(fieldname, "/NsIncomp2DFieldOutL%iC%iRe%.0fmt%.0f.mat",L,chi,Re,round(t*1000));
            strcat(currentSavename, fieldname );
            tntNetworksSave(currentSavename, u1, u2);
        }
    }
    
    /* Save Kinetic energy and Enstrophy as functions of time */
    currentSavename[0] = '\0';
    strcat(currentSavename,savepath);
    sprintf(txtfilename, "/NsIncomp2DEkinEnstL%iC%iRe%.0f.mat",L,chi,Re);
    strcat(currentSavename, txtfilename);
    printf("\nSaving kinetic energy and enstrophy at t=%f.",T); fflush(stdout);
    tntDoubleArraysSave(currentSavename, Ekin, Enst);
    tntDoubleParamsSave(currentSavename, Re,T,dt, lengthScale);
    
    /* Free everything and exit */
    tntDoubleArrayFree(&Ekin);
    tntDoubleArrayFree(&Enst);

    tntNetworkFree(&u1);
    tntNetworkFree(&u2);
    
    tntNetworkFree(&d1);
    tntNetworkFree(&d2);
    
    tntNetworkFree(&dd1);
    tntNetworkFree(&dd2);
    
    tntNetworkFree(&d1xd1);
    tntNetworkFree(&d1xd2);
    tntNetworkFree(&d2xd2);
    
    tntNodeArrayFree(&R1d1xd1R1_left); tntNodeArrayFree(&R1d1xd1R1_right);
    tntNodeArrayFree(&R1d1xd2R2_left); tntNodeArrayFree(&R1d1xd2R2_right);
    tntNodeArrayFree(&R2d2xd2R2_left); tntNodeArrayFree(&R2d2xd2R2_right);
    
}

