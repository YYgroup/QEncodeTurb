/*
 Authors: Nikita Gourianov
 (c) Nikita Gourianov 2019
*/

/* Include header for NS and PDE functions */
#include "tntPde.h"
#include "tntIncompNs.h"

/*! \ingroup ns
 * Performs a single, optimising sweep of the MPS (explicit time-stepping). *
 * \return No return value.
 */

void tntMpsNsIncomp3D1sSweep(tntNetwork v1, /*!< Assumes canonical centre at site "start_site" */
                             tntNetwork v2, /*!< Assumes canonical centre at site "start_site" */
                             tntNetwork v3, /*!< Assumes canonical centre at site "start_site" */
                             tntNetwork u1, tntNetwork u2, tntNetwork u3, tntNetwork uMidpoint1, tntNetwork uMidpoint2, tntNetwork uMidpoint3, tntNetwork* BurgMpoSum, 
                             tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d1xd3,
                             tntNetwork d2xd2, tntNetwork d2xd3, tntNetwork d3xd3,
                             tntNodeArray *R1xu1_left, tntNodeArray *R1xu1_right, 
                             tntNodeArray **R1xBxumid1_leftSum, tntNodeArray **R1xBxumid1_rightSum, 
                             tntNodeArray *R2xu2_left, tntNodeArray *R2xu2_right, 
                             tntNodeArray **R2xBxumid2_leftSum, tntNodeArray **R2xBxumid2_rightSum, 
                             tntNodeArray *R3xu3_left, tntNodeArray *R3xu3_right, 
                             tntNodeArray **R3xBxumid3_leftSum, tntNodeArray **R3xBxumid3_rightSum, 
                             tntNodeArray *R1d1xd1R1_left, tntNodeArray *R1d1xd1R1_right, 
                             tntNodeArray *R1d1xd2R2_left, tntNodeArray *R1d1xd2R2_right, 
                             tntNodeArray *R1d1xd3R3_left, tntNodeArray *R1d1xd3R3_right, 
                             tntNodeArray *R2d2xd2R2_left, tntNodeArray *R2d2xd2R2_right, 
                             tntNodeArray *R2d2xd3R3_left, tntNodeArray *R2d2xd3R3_right, 
                             tntNodeArray *R3d3xd3R3_left, tntNodeArray *R3d3xd3R3_right, 
                             tntLegLabel dir, 
                             unsigned numBurgTerms, 
                             int start_site, 
                             int end_site, 
                             char preconditioner[]){
    
    
    /* Declare relevant variables */
    int L, j; /* Number of sites and site index. Here j is chosen to always refer to the left-most node in the pair of nodes being updated. */
    tntNode Q, R; //Needed to shift canonical centre
    tntNode V1, V2, V3; //Nodes to be time-evolved
    tntNode U1, U2, U3, UM1, UM2, UM3;// Used to update L-1 node
    tntNode dummy; //Dummy variable to be used when moving orthogonal centre
    char preconditioner_local[20]; //Specifies which preconditioning is actually to be done at the local node level.
    
    /* These variables are needed to send in the sum of advection-diffusion terms into the local tntNode update function tntIncomp3D1sLocalTE */
    tntNode *Bsum;
    tntNode *R1xBxumid1_leftSumValsAtj, *R1xBxumid1_rightSumValsAtj;
    tntNode *R2xBxumid2_leftSumValsAtj, *R2xBxumid2_rightSumValsAtj;
    tntNode *R3xBxumid3_leftSumValsAtj, *R3xBxumid3_rightSumValsAtj;
    tntNode D1xD1, D1xD2, D1xD3, D2xD2, D2xD3, D3xD3;
    
    L = tntMpsLength(v1);
    
    /* Allocate memory to the node arrays */
    Bsum = malloc(numBurgTerms*sizeof(tntNode));
    R1xBxumid1_leftSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R1xBxumid1_rightSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R2xBxumid2_leftSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R2xBxumid2_rightSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R3xBxumid3_leftSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    R3xBxumid3_rightSumValsAtj = malloc(numBurgTerms*sizeof(tntNode));
    
    if ('R' == dir[0]) {
        /* Set local preconditioner to be Right block-inverse */
        if( strcmp(preconditioner,"None") ){
            sprintf(preconditioner_local,"Right Block-Inverse");
        }
        else{
            sprintf(preconditioner_local,"None");
        }
        
        /* Set the current tensors to be the first of their respective networks, except the advection-diffusion terms */
        V1 = tntNodeFindFirst(v1);
        V2 = tntNodeFindFirst(v2);
        V3 = tntNodeFindFirst(v3);
        U1 = tntNodeFindFirst(u1);
        U2 = tntNodeFindFirst(u2);
        U3 = tntNodeFindFirst(u3);
        UM1 = tntNodeFindFirst(uMidpoint1);
        UM2 = tntNodeFindFirst(uMidpoint2);
        UM3 = tntNodeFindFirst(uMidpoint3);
        D1xD1 = tntNodeFindFirst(d1xd1);
        D1xD2 = tntNodeFindFirst(d1xd2);
        D1xD3 = tntNodeFindFirst(d1xd3);
        D2xD2 = tntNodeFindFirst(d2xd2);
        D2xD3 = tntNodeFindFirst(d2xd3);
        D3xD3 = tntNodeFindFirst(d3xd3);
        
        /* Move forward until current tensors are start tensors, except the advection-diffusion terms */
        for (j = 0; j < start_site; j++) {
            V1 = tntNodeFindConn(V1,"R");
            V2 = tntNodeFindConn(V2,"R");
            V3 = tntNodeFindConn(V3,"R");
            U1 = tntNodeFindConn(U1,"R");
            U2 = tntNodeFindConn(U2,"R");
            U3 = tntNodeFindConn(U3,"R");
            UM1 = tntNodeFindConn(UM1,"R");
            UM2 = tntNodeFindConn(UM2,"R");
            UM3 = tntNodeFindConn(UM3,"R");
            
            D1xD1 = tntNodeFindConn(D1xD1,"R");
            D1xD2 = tntNodeFindConn(D1xD2,"R");
            D1xD3 = tntNodeFindConn(D1xD3,"R");
            D2xD2 = tntNodeFindConn(D2xD2,"R");
            D2xD3 = tntNodeFindConn(D2xD3,"R");
            D3xD3 = tntNodeFindConn(D3xD3,"R");
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
                R3xBxumid3_leftSumValsAtj[i] = ((*R3xBxumid3_leftSum)[i]).vals[j];
                R3xBxumid3_rightSumValsAtj[i] = ((*R3xBxumid3_rightSum)[i]).vals[j];
            }
            
            /* Then perform a 1-site local update */
            tntNodeNsIncomp3D1sLocalUpdate(&V1, &V2, &V3,
                                           U1, U2, U3, UM1, UM2, UM3, Bsum, D1xD1, D1xD2, D1xD3, D2xD2, D2xD3, D3xD3,
                                           R1xu1_left->vals[j], R1xu1_right->vals[j], R1xBxumid1_leftSumValsAtj, R1xBxumid1_rightSumValsAtj,
                                           R1d1xd1R1_left->vals[j], R1d1xd1R1_right->vals[j], R1d1xd2R2_left->vals[j], R1d1xd2R2_right->vals[j],
                                           R1d1xd3R3_left->vals[j], R1d1xd3R3_right->vals[j],
                                           R2xu2_left->vals[j], R2xu2_right->vals[j], R2xBxumid2_leftSumValsAtj, R2xBxumid2_rightSumValsAtj,
                                           R2d2xd2R2_left->vals[j], R2d2xd2R2_right->vals[j], R2d2xd3R3_left->vals[j], R2d2xd3R3_right->vals[j],
                                           R3xu3_left->vals[j], R3xu3_right->vals[j], R3xBxumid3_leftSumValsAtj, R3xBxumid3_rightSumValsAtj,
                                           R3d3xd3R3_left->vals[j], R3d3xd3R3_right->vals[j],
                                           numBurgTerms, preconditioner_local);
            
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
            
            /* Shift canonical centre from j to j+1 of field 3 */
            Q = tntNodeQR(V3,"LD", "R", "L");
            R = tntNodeFindConn(Q,"R");
            V3 = tntNodeContract(R, tntNodeFindConn(R,"R")); //Canonical centre now at site j+1
            V3 = tntNodeFindConn(V3,"L"); //Go back to site j
            
            /* Update node arrays */
            /* Start with updating the MPO ones */
            tntMpsPdeMpoUpdate(V1, D1xD1, V1, "R", j, R1d1xd1R1_left);
            tntMpsPdeMpoUpdate(V2, D1xD2, V1, "R", j, R1d1xd2R2_left);
            tntMpsPdeMpoUpdate(V3, D1xD3, V1, "R", j, R1d1xd3R3_left);
            tntMpsPdeMpoUpdate(V2, D2xD2, V2, "R", j, R2d2xd2R2_left);
            tntMpsPdeMpoUpdate(V3, D2xD3, V2, "R", j, R2d2xd3R3_left);
            tntMpsPdeMpoUpdate(V3, D3xD3, V3, "R", j, R3d3xd3R3_left);

            /* Then update the MPS ones between v and u fields (between updated and original fields) */
            tntMpsPdeMpsUpdate(U1, V1, "R", j, R1xu1_left);
            tntMpsPdeMpsUpdate(U2, V2, "R", j, R2xu2_left);
            tntMpsPdeMpsUpdate(U3, V3, "R", j, R3xu3_left);
            
            /* Finally, update the MPS ones between v and the advection-diffusion terms */
            for(unsigned i = 0; i<numBurgTerms; i++){
                tntMpsPdeMpoUpdate(UM1, Bsum[i], V1, "R", j, &((*R1xBxumid1_leftSum)[i]) );
                tntMpsPdeMpoUpdate(UM2, Bsum[i], V2, "R", j, &((*R2xBxumid2_leftSum)[i]) );
                tntMpsPdeMpoUpdate(UM3, Bsum[i], V3, "R", j, &((*R3xBxumid3_leftSum)[i]) );
            }
            
            /* Go to next set of nodes for all the MPSes*/
            V1 = tntNodeFindConn(V1,"R");
            V2 = tntNodeFindConn(V2,"R");
            V3 = tntNodeFindConn(V3,"R");
            U1 = tntNodeFindConn(U1,"R");
            U2 = tntNodeFindConn(U2,"R");
            U3 = tntNodeFindConn(U3,"R");
            UM1 = tntNodeFindConn(UM1,"R");
            UM2 = tntNodeFindConn(UM2,"R");
            UM3 = tntNodeFindConn(UM3,"R");
            D1xD1 = tntNodeFindConn(D1xD1,"R");
            D1xD2 = tntNodeFindConn(D1xD2,"R");
            D1xD3 = tntNodeFindConn(D1xD3,"R");
            D2xD2 = tntNodeFindConn(D2xD2,"R");
            D2xD3 = tntNodeFindConn(D2xD3,"R");
            D3xD3 = tntNodeFindConn(D3xD3,"R");
            for(unsigned i = 0; i<numBurgTerms; i++){
                Bsum[i] = tntNodeFindConn(Bsum[i],"R");
            }
            
        }
    }
    else{
        /* Set local preconditioner to be Left block-inverse */
        if( strcmp(preconditioner,"None") ){
            sprintf(preconditioner_local,"Left Block-Inverse");
        }
        else{
            sprintf(preconditioner_local,"None");
        }
        
        /* Set the current tensors to be the last of their respective networks, except the advection-diffusion terms */
        V1 = tntNodeFindLast(v1);
        V2 = tntNodeFindLast(v2);
        V3 = tntNodeFindLast(v3);
        U1 = tntNodeFindLast(u1);
        U2 = tntNodeFindLast(u2);
        U3 = tntNodeFindLast(u3);
        UM1 = tntNodeFindLast(uMidpoint1);
        UM2 = tntNodeFindLast(uMidpoint2);
        UM3 = tntNodeFindLast(uMidpoint3);
        D1xD1 = tntNodeFindLast(d1xd1);
        D1xD2 = tntNodeFindLast(d1xd2);
        D1xD3 = tntNodeFindLast(d1xd3);
        D2xD2 = tntNodeFindLast(d2xd2);
        D2xD3 = tntNodeFindLast(d2xd3);
        D3xD3 = tntNodeFindLast(d3xd3);

        /* Move backwards until current tensors are start tensors, except the advection-diffusion terms */
        for (j = L-1; j > start_site; j--) {
            V1 = tntNodeFindConn(V1,"L");
            V2 = tntNodeFindConn(V2,"L");
            V3 = tntNodeFindConn(V3,"L");
            U1 = tntNodeFindConn(U1,"L");
            U2 = tntNodeFindConn(U2,"L");
            U3 = tntNodeFindConn(U3,"L");
            UM1 = tntNodeFindConn(UM1,"L");
            UM2 = tntNodeFindConn(UM2,"L");
            UM3 = tntNodeFindConn(UM3,"L");
            D1xD1 = tntNodeFindConn(D1xD1,"L");
            D1xD2 = tntNodeFindConn(D1xD2,"L");
            D1xD3 = tntNodeFindConn(D1xD3,"L");
            D2xD2 = tntNodeFindConn(D2xD2,"L");
            D2xD3 = tntNodeFindConn(D2xD3,"L");
            D3xD3 = tntNodeFindConn(D3xD3,"L");
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
                R3xBxumid3_leftSumValsAtj[i] = ((*R3xBxumid3_leftSum)[i]).vals[j];
                R3xBxumid3_rightSumValsAtj[i] = ((*R3xBxumid3_rightSum)[i]).vals[j];
                
            }
  
            
            tntNodeNsIncomp3D1sLocalUpdate(&V1, &V2, &V3,
                                           U1, U2, U3, UM1, UM2, UM3, Bsum, D1xD1, D1xD2, D1xD3, D2xD2, D2xD3, D3xD3,
                                           R1xu1_left->vals[j], R1xu1_right->vals[j], R1xBxumid1_leftSumValsAtj, R1xBxumid1_rightSumValsAtj,
                                           R1d1xd1R1_left->vals[j], R1d1xd1R1_right->vals[j], R1d1xd2R2_left->vals[j], R1d1xd2R2_right->vals[j],
                                           R1d1xd3R3_left->vals[j], R1d1xd3R3_right->vals[j],
                                           R2xu2_left->vals[j], R2xu2_right->vals[j], R2xBxumid2_leftSumValsAtj, R2xBxumid2_rightSumValsAtj,
                                           R2d2xd2R2_left->vals[j], R2d2xd2R2_right->vals[j], R2d2xd3R3_left->vals[j], R2d2xd3R3_right->vals[j],
                                           R3xu3_left->vals[j], R3xu3_right->vals[j], R3xBxumid3_leftSumValsAtj, R3xBxumid3_rightSumValsAtj,
                                           R3d3xd3R3_left->vals[j], R3d3xd3R3_right->vals[j],
                                           numBurgTerms, preconditioner_local);
            
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
            
            /* Shift canonical centre from j to j-1 of field 2 */
            Q = tntNodeQR(V3,"RD", "L", "R");
            R = tntNodeFindConn(Q,"L");
            V3 = tntNodeContract(R, tntNodeFindConn(R,"L")); //Canonical centre now at site j-1
            V3 = tntNodeFindConn(V3,"R"); //Go back to site j
            
            /* Update node arrays */
            /* Start with updating the MPO ones */
            tntMpsPdeMpoUpdate(V1, D1xD1, V1, "L", j, R1d1xd1R1_right);
            tntMpsPdeMpoUpdate(V2, D1xD2, V1, "L", j, R1d1xd2R2_right);
            tntMpsPdeMpoUpdate(V3, D1xD3, V1, "L", j, R1d1xd3R3_right);
            tntMpsPdeMpoUpdate(V2, D2xD2, V2, "L", j, R2d2xd2R2_right);
            tntMpsPdeMpoUpdate(V3, D2xD3, V2, "L", j, R2d2xd3R3_right);
            tntMpsPdeMpoUpdate(V3, D3xD3, V3, "L", j, R3d3xd3R3_right);

            /* Then update the MPS ones between v and u fields (between updated and original fields) */
            tntMpsPdeMpsUpdate(U1,V1, "L", j, R1xu1_right);
            tntMpsPdeMpsUpdate(U2,V2, "L", j, R2xu2_right);
            tntMpsPdeMpsUpdate(U3,V3, "L", j, R3xu3_right);

            /* Finally, upadte the MPS ones between v and the advection-diffusion terms */
            for(unsigned i = 0; i<numBurgTerms; i++){
                tntMpsPdeMpoUpdate(UM1, Bsum[i],V1, "L", j, &((*R1xBxumid1_rightSum)[i]));
                tntMpsPdeMpoUpdate(UM2, Bsum[i],V2, "L", j, &((*R2xBxumid2_rightSum)[i]));
                tntMpsPdeMpoUpdate(UM3, Bsum[i],V3, "L", j, &((*R3xBxumid3_rightSum)[i]));
            }
            
            /* Go to next set of nodes for all the MPSes*/
            V1 = tntNodeFindConn(V1,"L");
            V2 = tntNodeFindConn(V2,"L");
            V3 = tntNodeFindConn(V3,"L");
            U1 = tntNodeFindConn(U1,"L");
            U2 = tntNodeFindConn(U2,"L");
            U3 = tntNodeFindConn(U3,"L");
            UM1 = tntNodeFindConn(UM1,"L");
            UM2 = tntNodeFindConn(UM2,"L");
            UM3 = tntNodeFindConn(UM3,"L");
            D1xD1 = tntNodeFindConn(D1xD1,"L");
            D1xD2 = tntNodeFindConn(D1xD2,"L");
            D1xD3 = tntNodeFindConn(D1xD3,"L");
            D2xD2 = tntNodeFindConn(D2xD2,"L");
            D2xD3 = tntNodeFindConn(D2xD3,"L");
            D3xD3 = tntNodeFindConn(D3xD3,"L");
            
            for(unsigned i = 0; i<numBurgTerms; i++){
                Bsum[i] = tntNodeFindConn(Bsum[i],"L");
            }
        }
    }
    
    /* Must time-evolve final site as well. */
    
    /* Set local preconditioner */
    if( strcmp(preconditioner,"None") ){
        if(end_site == 0){
            sprintf(preconditioner_local,"Right Block-Inverse"); }
        else if(end_site == L-1){
            sprintf(preconditioner_local,"Left Block-Inverse"); }
        else{
            sprintf(preconditioner_local,"None");
            printf("\nWarning: this line shouldn't trigger in tntMps1sNsIncompImex3D. Investigate.\n"); }
    }
    else{
        sprintf(preconditioner_local,"None");
    }
    
    /* Gather the precontracted advection-diffusion networks into arrays before time-evolving */
    for(unsigned i = 0; i<numBurgTerms; i++){
        R1xBxumid1_leftSumValsAtj[i] = ((*R1xBxumid1_leftSum)[i]).vals[end_site];
        R1xBxumid1_rightSumValsAtj[i] = ((*R1xBxumid1_rightSum)[i]).vals[end_site];
        R2xBxumid2_leftSumValsAtj[i] = ((*R2xBxumid2_leftSum)[i]).vals[end_site];
        R2xBxumid2_rightSumValsAtj[i] = ((*R2xBxumid2_rightSum)[i]).vals[end_site];
        R3xBxumid3_leftSumValsAtj[i] = ((*R3xBxumid3_leftSum)[i]).vals[end_site];
        R3xBxumid3_rightSumValsAtj[i] = ((*R3xBxumid3_rightSum)[i]).vals[end_site];
    }
    
    
    tntNodeNsIncomp3D1sLocalUpdate(&V1, &V2, &V3,
                               U1, U2, U3, UM1, UM2, UM3, Bsum, D1xD1, D1xD2, D1xD3, D2xD2, D2xD3, D3xD3,
                               R1xu1_left->vals[end_site], R1xu1_right->vals[end_site], R1xBxumid1_leftSumValsAtj, R1xBxumid1_rightSumValsAtj,
                               R1d1xd1R1_left->vals[end_site], R1d1xd1R1_right->vals[end_site], R1d1xd2R2_left->vals[end_site], R1d1xd2R2_right->vals[end_site],
                               R1d1xd3R3_left->vals[end_site], R1d1xd3R3_right->vals[end_site],
                               R2xu2_left->vals[end_site], R2xu2_right->vals[end_site], R2xBxumid2_leftSumValsAtj, R2xBxumid2_rightSumValsAtj,
                               R2d2xd2R2_left->vals[end_site], R2d2xd2R2_right->vals[end_site], R2d2xd3R3_left->vals[end_site], R2d2xd3R3_right->vals[end_site],
                               R3xu3_left->vals[end_site], R3xu3_right->vals[end_site], R3xBxumid3_leftSumValsAtj, R3xBxumid3_rightSumValsAtj,
                               R3d3xd3R3_left->vals[end_site], R3d3xd3R3_right->vals[end_site],
                               numBurgTerms, preconditioner_local);
 
    
    /* Free the arrays of nodes */
    free(Bsum);
    free(R1xBxumid1_leftSumValsAtj);
    free(R1xBxumid1_rightSumValsAtj);
    free(R2xBxumid2_leftSumValsAtj);
    free(R2xBxumid2_rightSumValsAtj);
    free(R3xBxumid3_leftSumValsAtj);
    free(R3xBxumid3_rightSumValsAtj);
}



/*! \ingroup ns
 * Performs a double-sweep optimising the MPS (explicit time-stepping).*
 * \return No return value.
 */

void tntMpsNsIncomp3D1sFullSweep(tntNetwork v1, tntNetwork v2, tntNetwork v3, 
                                 tntNetwork u1, tntNetwork u2, tntNetwork u3,  tntNetwork uMidpoint1, tntNetwork uMidpoint2, tntNetwork uMidpoint3, 
                                 tntNetwork *BurgMpoSum,  
                                 tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d1xd3, tntNetwork d2xd2, tntNetwork d2xd3, tntNetwork d3xd3,
                                 tntNodeArray *R1xu1_left, tntNodeArray *R1xu1_right, 
                                 tntNodeArray **R1xBxumid1_leftSum, tntNodeArray **R1xBxumid1_rightSum, 
                                 tntNodeArray *R2xu2_left, tntNodeArray *R2xu2_right, 
                                 tntNodeArray **R2xBxumid2_leftSum, tntNodeArray **R2xBxumid2_rightSum, 
                                 tntNodeArray *R3xu3_left, tntNodeArray *R3xu3_right, 
                                 tntNodeArray **R3xBxumid3_leftSum, tntNodeArray **R3xBxumid3_rightSum, 
                                 tntNodeArray *R1d1xd1R1_left, tntNodeArray *R1d1xd1R1_right,  
                                 tntNodeArray *R1d1xd2R2_left, tntNodeArray *R1d1xd2R2_right,   
                                 tntNodeArray *R1d1xd3R3_left, tntNodeArray *R1d1xd3R3_right,   
                                 tntNodeArray *R2d2xd2R2_left, tntNodeArray *R2d2xd2R2_right,  
                                 tntNodeArray *R2d2xd3R3_left, tntNodeArray *R2d2xd3R3_right,  
                                 tntNodeArray *R3d3xd3R3_left, tntNodeArray *R3d3xd3R3_right,  
                                 unsigned numBurgTerms, /* Number of terms in the advection-diffusion term 1 & 2 sums */
                                 char preconditioner[]){
    
    unsigned L = tntMpsLength(v1);
    
    if(L == 1){
        tntErrorPrint("Cannot integrate the incompressible NS equation forward if it just has one site! Exiting."); /* NO_COVERAGE */
        exit(-1);
    }
    else{
        /* Initialise relevant nodes */
        tntNode Q, R, V1, V2, V3; /* Needed to shift canonical centre from L-1 to L-2 */
        tntNode D1xD1, D1xD2, D1xD3, D2xD2, D2xD3, D3xD3;

        /* Sweep from node 0 to node L-1, i.e. from the very first to the very last */
        tntMpsNsIncomp3D1sSweep(v1, v2, v3,
                                u1, u2, u3,uMidpoint1, uMidpoint2, uMidpoint3, BurgMpoSum,
                                d1xd1, d1xd2, d1xd3, d2xd2, d2xd3, d3xd3,
                                R1xu1_left, R1xu1_right, R1xBxumid1_leftSum, R1xBxumid1_rightSum,
                                R2xu2_left, R2xu2_right, R2xBxumid2_leftSum, R2xBxumid2_rightSum,
                                R3xu3_left, R3xu3_right, R3xBxumid3_leftSum, R3xBxumid3_rightSum,
                                R1d1xd1R1_left, R1d1xd1R1_right, R1d1xd2R2_left, R1d1xd2R2_right, R1d1xd3R3_left, R1d1xd3R3_right,
                                R2d2xd2R2_left, R2d2xd2R2_right, R2d2xd3R3_left, R2d2xd3R3_right,
                                R3d3xd3R3_left, R3d3xd3R3_right,
                                "R", numBurgTerms, 0, L-1, preconditioner);

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

        /* Shift canonical centre from L-1 to L-2 of field 3 */
        V3 = tntNodeFindLast(v3);
        Q = tntNodeQR(V3,"RD", "L", "R");
        R = tntNodeFindConn(Q,"L");
        V3 = tntNodeContract(R, tntNodeFindConn(R,"L")); //Canonical centre now at site L-2
                
        /* Update node arrays */
        /* First define the required nodes at L-1 */
        V1 = tntNodeFindConn(V1,"R");
        V2 = tntNodeFindConn(V2,"R");
        V3 = tntNodeFindConn(V3,"R");
        
        D1xD1 = tntNodeFindLast(d1xd1);
        D1xD2 = tntNodeFindLast(d1xd2);
        D1xD3 = tntNodeFindLast(d1xd3);
        D2xD2 = tntNodeFindLast(d2xd2);
        D2xD3 = tntNodeFindLast(d2xd3);
        D3xD3 = tntNodeFindLast(d3xd3);
        
        /* Then do it for the sum of the advection-diffusion terms */
        for(unsigned i = 0; i<numBurgTerms; i++){
            tntMpsPdeMpoUpdate(tntNodeFindLast(uMidpoint1),tntNodeFindLast(BurgMpoSum[i]),V1, "L", L-1, &((*R1xBxumid1_rightSum)[i]) );
            tntMpsPdeMpoUpdate(tntNodeFindLast(uMidpoint2),tntNodeFindLast(BurgMpoSum[i]),V2, "L", L-1, &((*R2xBxumid2_rightSum)[i]) );
            tntMpsPdeMpoUpdate(tntNodeFindLast(uMidpoint3),tntNodeFindLast(BurgMpoSum[i]),V3, "L", L-1, &((*R3xBxumid3_rightSum)[i]) );
        }

        /* Then for the differencing networks: */
        tntMpsPdeMpoUpdate(V1, D1xD1, V1, "L", L-1, R1d1xd1R1_right);
        tntMpsPdeMpoUpdate(V2, D1xD2, V1, "L", L-1, R1d1xd2R2_right);
        tntMpsPdeMpoUpdate(V3, D1xD3, V1, "L", L-1, R1d1xd3R3_right);
        tntMpsPdeMpoUpdate(V2, D2xD2, V2, "L", L-1, R2d2xd2R2_right);
        tntMpsPdeMpoUpdate(V3, D2xD3, V2, "L", L-1, R2d2xd3R3_right);
        tntMpsPdeMpoUpdate(V3, D3xD3, V3, "L", L-1, R3d3xd3R3_right);


        /* And finally for the initial fields: */
        tntMpsPdeMpsUpdate(tntNodeFindLast(u1),V1, "L", L-1, R1xu1_right);
        tntMpsPdeMpsUpdate(tntNodeFindLast(u2),V2, "L", L-1, R2xu2_right);
        tntMpsPdeMpsUpdate(tntNodeFindLast(u3),V3, "L", L-1, R3xu3_right);
        
        /* Now Sweep from node L-2 to node 0, i.e. from the second-to-last to the first, and finish */
        tntMpsNsIncomp3D1sSweep(v1, v2, v3,
                                u1, u2, u3, uMidpoint1, uMidpoint2, uMidpoint3, BurgMpoSum,
                                d1xd1, d1xd2, d1xd3, d2xd2, d2xd3, d3xd3,
                                R1xu1_left, R1xu1_right, R1xBxumid1_leftSum, R1xBxumid1_rightSum,
                                R2xu2_left, R2xu2_right, R2xBxumid2_leftSum, R2xBxumid2_rightSum,
                                R3xu3_left, R3xu3_right, R3xBxumid3_leftSum, R3xBxumid3_rightSum,
                                R1d1xd1R1_left, R1d1xd1R1_right, R1d1xd2R2_left, R1d1xd2R2_right, R1d1xd3R3_left, R1d1xd3R3_right,
                                R2d2xd2R2_left, R2d2xd2R2_right, R2d2xd3R3_left, R2d2xd3R3_right,
                                R3d3xd3R3_left, R3d3xd3R3_right,
                                "L", numBurgTerms, L-2, 0, preconditioner);
    }
    
}

/*! \ingroup ns
 * Performs full Runge-Kutta 2 time-step  of the MPS (explicit time-stepping).
 * \return No return value.
 */

void tntMpsNsIncomp3D1sTimeStep(tntNetwork v1, tntNetwork v2, tntNetwork v3,  
                            tntNetwork u1, tntNetwork u2, tntNetwork u3, tntNetwork uMidpoint1, tntNetwork uMidpoint2, tntNetwork uMidpoint3, 
                            tntNetwork *BurgMpoSum,
                            tntNetwork d1xd1, tntNetwork d1xd2, tntNetwork d1xd3, tntNetwork d2xd2, tntNetwork d2xd3, tntNetwork d3xd3, 
                            tntNodeArray *R1d1xd1R1_left, tntNodeArray *R1d1xd1R1_right,  
                            tntNodeArray *R1d1xd2R2_left, tntNodeArray *R1d1xd2R2_right,   
                            tntNodeArray *R1d1xd3R3_left, tntNodeArray *R1d1xd3R3_right,  
                            tntNodeArray *R2d2xd2R2_left, tntNodeArray *R2d2xd2R2_right,  
                            tntNodeArray *R2d2xd3R3_left, tntNodeArray *R2d2xd3R3_right,   
                            tntNodeArray *R3d3xd3R3_left, tntNodeArray *R3d3xd3R3_right,  
                            unsigned numBurgTerms, /* Number of terms in the Burger MPOs sum array */
                            double epsilon, /* convergence number for sweeping */
                            char dynamicAdjustments[]){
    
    /* Declare relevant variables */
    int numRuns;
    double Ek0, Ek1 = 0; //Kinetic energy of field
    char preconditioner[20];
    unsigned L = tntMpsLength(v1);
    tntNodeArray R1xu1_left, R1xu1_right;
    tntNodeArray R2xu2_left, R2xu2_right;
    tntNodeArray R3xu3_left, R3xu3_right;
    tntNodeArray *R1xBxumid1_leftSum, *R1xBxumid1_rightSum;
    tntNodeArray *R2xBxumid2_leftSum, *R2xBxumid2_rightSum;
    tntNodeArray *R3xBxumid3_leftSum, *R3xBxumid3_rightSum;
    tntNetwork burgMpo;

    /* Allocate memory to the arrays of node arrays */
    R1xBxumid1_leftSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R1xBxumid1_rightSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R2xBxumid2_leftSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R2xBxumid2_rightSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R3xBxumid3_leftSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    R3xBxumid3_rightSum = malloc(numBurgTerms*sizeof(tntNodeArray));
    
    /* Produce the pre-contracted networks */
    tntMpsPdeMpsInit(u1,v1,&R1xu1_left,&R1xu1_right);
    tntMpsPdeMpsInit(u2,v2,&R2xu2_left,&R2xu2_right);
    tntMpsPdeMpsInit(u3,v3,&R3xu3_left,&R3xu3_right);
    
    for(unsigned i = 0; i<numBurgTerms; i++){
        burgMpo = BurgMpoSum[i];
        tntMpsPdeMpoInit(uMidpoint1,burgMpo,v1,&(R1xBxumid1_leftSum[i]),&(R1xBxumid1_rightSum[i]));
        tntMpsPdeMpoInit(uMidpoint2,burgMpo,v2,&(R2xBxumid2_leftSum[i]),&(R2xBxumid2_rightSum[i]));
        tntMpsPdeMpoInit(uMidpoint3,burgMpo,v3,&(R3xBxumid3_leftSum[i]),&(R3xBxumid3_rightSum[i]));
    }
    
    /* Sweep until convergence */
    numRuns = 0;
    Ek0 = 0;
    Ek1 = 2*epsilon;
    while( fabs((Ek1-Ek0)/Ek0) > epsilon || numRuns<2){
        
        /* Preconditioner is turned off by default */
        sprintf(preconditioner,"None");
        /* Adjust simulation settings dynamically based on numRuns, if requested */
        if( strcmp(dynamicAdjustments,"None") ){
            /* If taking too long, turn on preconditioner */
            if( ((numRuns+1)%500) == 0 ){ //Precon doesnt work properly, as it leaks memory.. So disabling it for now.
                sprintf(preconditioner,"Yes");
                printf("\nTaking too long. Turning on preconditioner in an attempt to speed things up.");
            }
            /* Give up after 50 runs */
            else if(numRuns>50){
                printf("\nBeen over a 50 runs now. It's probably not going to converge, so giving up.");
                break;
            }
        }
        
        tntMpsNsIncomp3D1sFullSweep(v1, v2, v3,
                                    u1, u2, u3, uMidpoint1, uMidpoint2, uMidpoint3,
                                    BurgMpoSum,
                                    d1xd1, d1xd2, d1xd3, d2xd2, d2xd3, d3xd3,
                                    &R1xu1_left, &R1xu1_right,
                                    &R1xBxumid1_leftSum, &R1xBxumid1_rightSum,
                                    &R2xu2_left, &R2xu2_right,
                                    &R2xBxumid2_leftSum, &R2xBxumid2_rightSum,
                                    &R3xu3_left, &R3xu3_right,
                                    &R3xBxumid3_leftSum, &R3xBxumid3_rightSum,
                                    R1d1xd1R1_left, R1d1xd1R1_right,
                                    R1d1xd2R2_left, R1d1xd2R2_right,
                                    R1d1xd3R3_left, R1d1xd3R3_right,
                                    R2d2xd2R2_left, R2d2xd2R2_right,
                                    R2d2xd3R3_left, R2d2xd3R3_right,
                                    R3d3xd3R3_left, R3d3xd3R3_right,
                                    numBurgTerms,preconditioner);
        
        numRuns++;
        
        /* Update Kinetic energy */
        Ek0=Ek1;
        Ek1 = tntMpsSelfProduct(v1,-1)+tntMpsSelfProduct(v2,-1)+tntMpsSelfProduct(v3,-1);
        printf("\nDouble-Sweep completed. Ek1=%.6e, (Ek1-Ek0)/Ek0=%.6e", Ek1, (Ek1-Ek0)/Ek0 );
    }
    
    printf("\n=====Convergence reached at epsilon=%.5e after %i runs=====",epsilon, numRuns);
    
    /* Free the node arrays */
    tntNodeArrayFree(&R1xu1_left); tntNodeArrayFree(&R1xu1_right);
    tntNodeArrayFree(&R2xu2_left); tntNodeArrayFree(&R2xu2_right);
    tntNodeArrayFree(&R3xu3_left); tntNodeArrayFree(&R3xu3_right);
    
    for(unsigned i = 0; i<numBurgTerms; i++){
        tntNodeArrayFree( &(R1xBxumid1_leftSum[i]));
        tntNodeArrayFree( &(R1xBxumid1_rightSum[i]));
        tntNodeArrayFree( &(R2xBxumid2_leftSum[i]));
        tntNodeArrayFree( &(R2xBxumid2_rightSum[i]));
        tntNodeArrayFree( &(R3xBxumid3_leftSum[i]));
        tntNodeArrayFree( &(R3xBxumid3_rightSum[i]));
    }

    /* Free the array of node arrays */
    free(R1xBxumid1_leftSum); free(R1xBxumid1_rightSum);
    free(R2xBxumid2_leftSum); free(R2xBxumid2_rightSum);
    free(R3xBxumid3_leftSum); free(R3xBxumid3_rightSum);
    
}


/*! \ingroup ns
 * Performs multiple Runge-Kutta 2 time-steps of the MPS (explicit time-stepping).
 * \return No return value.
 */

void tntMpsBoxTurbulence3DTimeEvolutionRK2(tntNetwork *v1, tntNetwork *v2, tntNetwork *v3,  /* Inputed MPS flow field  */
                                           double t0, /* Starting time */
                                           double dt,  /* Time step */
                                           double T,  /* Final time */
                                           double T_sampling,    /* Sampling period (i.e. 1/T_sampling = sampling frequency) */
                                           char savepath[TNT_STRLEN], /* Output path */
                                           char BC1[], char BC2[], char BC3[], /* Boundary conditions of the box */
                                           double Re,    /* Reynolds number */
                                           double lengthScale, /* Length of box */
                                           double lambda,  /* Penalty coefficient */
                                           double epsilon,  /* Convergence sensitivity */
                                           int chi, int chi_mpo, /* Bond dimension for MPSes and MPOs, respectively */
                                           double *errp /* Dummy variable */ ){
    
    extern tntSystem systype;
    
    /* First of all, make sure the user has turned of symmetries. */
    if(TNT_SYMM_NONE != tntSymmTypeGet()){
        char warning_msg[TNT_STRLEN];
        strcpy(warning_msg,"tntMpsBoxTurbulence3DTimeEvolutionRK2 doesn't work with symmetries! Turning them OFF and continuing!");
        tntSymmTypeUnset();
        tntWarningPrint(warning_msg);
    }
    
    /* Declare relevant variables */
    char currentSavename[TNT_STRLEN];
    char fieldname[TNT_STRLEN], txtfilename[TNT_STRLEN], dynamicAdjustments[20];;
    tntDoubleArray Ekin, Enst;
    tntDoubleArray mpsEntanSpec1, mpsEntanSpec2, mpsEntanSpec3, mpsEntanEntr1, mpsEntanEntr2, mpsEntanEntr3;
    tntComplexArray dummy_array;
    tntNetwork u1, u2, u3, uMidpoint1, uMidpoint2, uMidpoint3;
    tntNetwork *BurgMpoSum;
    tntNetwork d1, d2, d3;
    tntNetwork dd1, dd2, dd3;
    tntNetwork d1xd1, d1xd2, d1xd3, d2xd2, d2xd3, d3xd3;
    tntNodeArray R1d1xd1R1_left, R1d1xd1R1_right;
    tntNodeArray R1d1xd2R2_left, R1d1xd2R2_right;
    tntNodeArray R1d1xd3R3_left, R1d1xd3R3_right;
    tntNodeArray R2d2xd2R2_left, R2d2xd2R2_right;
    tntNodeArray R2d2xd3R3_left, R2d2xd3R3_right;
    tntNodeArray R3d3xd3R3_left, R3d3xd3R3_right;
    tntNode A, B, C;
    double chi_Lleg, chi_Rleg;
    double t;
    double h, Ek0, Ek;
    unsigned numBurgTerms;
    unsigned L = tntMpsLength(*v1);
    int step, stepsBetweenSampling, L_dim, R_dim;
    
    /*  unless simulation is just continuing from a previous checkpoint (when t0!=0), contract the original field down to the appropriate bond dimension and use it as the "starting guess" for the time-evolved field. Also, put it in appropriate orthogonal form. */
    if(t0==0){
    	tntMpsOrth(*v1,L-1);
    	tntMpsOrth(*v2,L-1);
    	tntMpsOrth(*v3,L-1);
    	*errp = *errp + tntMpsTruncate(*v1, 0, chi);
    	*errp = *errp + tntMpsTruncate(*v2, 0, chi);
    	*errp = *errp + tntMpsTruncate(*v3, 0, chi);
    }
    /* Initialise differencing MPO & operators*/
    d1 = tntMpsPde3DCreateOrder1Diff1(BC1,L);
    d2 = tntMpsPde3DCreateOrder1Diff2(BC2,L);
    d3 = tntMpsPde3DCreateOrder1Diff3(BC3,L);
    dd1 = tntMpsPde3DCreateOrder2Diff1(BC1,L);
    dd2 = tntMpsPde3DCreateOrder2Diff2(BC2,L);
    dd3 = tntMpsPde3DCreateOrder2Diff3(BC3,L);
    
    /* Create MPO networks as products between differencers and their adjoints */
    d1xd1 = tntMpoMpoProduct(d1, d1);
    d1xd2 = tntMpoMpoProduct(d2, d1);
    d1xd3 = tntMpoMpoProduct(d3, d1);
    d2xd2 = tntMpoMpoProduct(d2, d2);
    d2xd3 = tntMpoMpoProduct(d3, d2);
    d3xd3 = tntMpoMpoProduct(d3, d3);

    /* Exploit sparseness to decrease bond dimension */
    *errp += tntMpoTruncate(d1xd1, -1);
    *errp += tntMpoTruncate(d1xd2, -1);
    *errp += tntMpoTruncate(d1xd3, -1);
    *errp += tntMpoTruncate(d2xd2, -1);
    *errp += tntMpoTruncate(d2xd3, -1);
    *errp += tntMpoTruncate(d3xd3, -1);
    
    /* Scale by -1 to account for the fact that the differencing operator is anti-symmetric, and by lambda. */
    tntNodeScaleReal(tntNodeFindFirst(d1xd1),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d1xd2),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d1xd3),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d2xd2),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d2xd3),-lambda);
    tntNodeScaleReal(tntNodeFindFirst(d3xd3),-lambda);
        
    /* Scale by lengthScale to adjust for domain size */
    tntNodeScaleReal(tntNodeFindFirst(d1),1/lengthScale);
    tntNodeScaleReal(tntNodeFindFirst(d2),1/lengthScale);
    tntNodeScaleReal(tntNodeFindFirst(d3),1/lengthScale);
    
    tntNodeScaleReal(tntNodeFindFirst(dd1),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(dd2),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(dd3),1/pow(lengthScale,2));
    
    tntNodeScaleReal(tntNodeFindFirst(d1xd1),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(d1xd2),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(d1xd3),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(d2xd2),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(d2xd3),1/pow(lengthScale,2));
    tntNodeScaleReal(tntNodeFindFirst(d3xd3),1/pow(lengthScale,2));
    
    /* Prepare system for time-evolution */
    t=t0; //Start at time = t0.
    step = 0;
    stepsBetweenSampling = floor(T_sampling/dt);
    if(stepsBetweenSampling == 0){
        stepsBetweenSampling = 1;
    }
    
    /* Initialise initial field */
    u1 = tntNetworkCopy(*v1);
    u2 = tntNetworkCopy(*v2);
    u3 = tntNetworkCopy(*v3);
    
    /* Set initial fields to be random to support prolongation, unless simulation is just continuing from a previous checkpoint (when t0!=0) */
    if(t0==0){
    	tntNetworkFree(v1);
   	 tntNetworkFree(v2);
   	 tntNetworkFree(v3);
	
	    *v1 = tntNetworkCreate();
	    *v2 = tntNetworkCreate();
	    *v3 = tntNetworkCreate();
	    R_dim = 1;
	    for(unsigned i = 0;i<L; i++){
	       //  Set dimensions
	        L_dim = R_dim;
       	 	R_dim = MIN(MIN(pow(8,i+1), pow(8,L-i-1)), chi);
	
        	// Initialise nodes to be zeros
        	dummy_array = tntComplexArrayAlloc(L_dim*8,R_dim);

        	A = tntNodeCreate(&dummy_array,"DLR",8,L_dim, R_dim);
       		 B = tntNodeCreate(&dummy_array,"DLR",8,L_dim, R_dim);
        	C = tntNodeCreate(&dummy_array,"DLR",8,L_dim, R_dim);

        	// Insert nodes
        	tntNodeInsertAtEnd(A,"L","R",*v1);
        	tntNodeInsertAtEnd(B,"L","R",*v2);
        	tntNodeInsertAtEnd(C,"L","R",*v3);
    	    }
        tntMpsOrth(*v1,0);
        tntMpsOrth(*v2,0);
        tntMpsOrth(*v3,0);
    }

    /* Initialise Kinetic energy and Enstrophy arrays */
    Ekin = tntDoubleArrayAlloc(ceil( (T-t0)/dt)+1);
    Enst = tntDoubleArrayAlloc(ceil( (T-t0)/dt)+1);
    
    /* Initialise pre-contracted networks */
    tntMpsPdeMpoInit(*v1,d1xd1,*v1,&R1d1xd1R1_left,&R1d1xd1R1_right);
    tntMpsPdeMpoInit(*v2,d1xd2,*v1,&R1d1xd2R2_left,&R1d1xd2R2_right);
    tntMpsPdeMpoInit(*v3,d1xd3,*v1,&R1d1xd3R3_left,&R1d1xd3R3_right);
    tntMpsPdeMpoInit(*v2,d2xd2,*v2,&R2d2xd2R2_left,&R2d2xd2R2_right);
    tntMpsPdeMpoInit(*v3,d2xd3,*v2,&R2d2xd3R3_left,&R2d2xd3R3_right);
    tntMpsPdeMpoInit(*v3,d3xd3,*v3,&R3d3xd3R3_left,&R3d3xd3R3_right);
    
    /* Initialise arrays containing the entanglement spectrum and entropy */
    mpsEntanSpec1 = tntDoubleArrayAlloc(chi, L-1); mpsEntanSpec2 = tntDoubleArrayAlloc(chi, L-1); mpsEntanSpec3 = tntDoubleArrayAlloc(chi, L-1);
    mpsEntanEntr1 = tntDoubleArrayAlloc(1, L-1); mpsEntanEntr2 = tntDoubleArrayAlloc(1, L-1); mpsEntanEntr3 = tntDoubleArrayAlloc(1, L-1);
    
    /* Calculate initial entanglement spectrum and entropy. */
    tntMpsCalcEntanSpec(u1, chi, 0, &mpsEntanSpec1, &mpsEntanEntr1);
    tntMpsCalcEntanSpec(u2, chi, 0, &mpsEntanSpec2, &mpsEntanEntr2);
    tntMpsCalcEntanSpec(u3, chi, 0, &mpsEntanSpec3, &mpsEntanEntr3);
    
    /* Save initial field */
    printf("\nSaving initial fields at t=%f.",0);
    printf("\nNormalised Entanglement Entropies at every bond = ");
    for(int i = 0; i<L-1; i++){
        printf("[%.4f,%.4f,%.4f]_%i; ",mpsEntanEntr1.vals[i],mpsEntanEntr2.vals[i],mpsEntanEntr3.vals[i],i+1);
    } printf("\n");     fflush(stdout);
    
    currentSavename[0] = '\0'; strcat(currentSavename,savepath);
    sprintf(fieldname, "/NsIncomp3DFieldOutL%iC%iRe%.0fmt%.0f.mat",L,chi,Re,round(t*1000));
    strcat(currentSavename, fieldname );
    tntNetworksSave(currentSavename, u1, u2, u3);
    tntDoubleParamsSave(currentSavename, Re,t,T_sampling, lengthScale);
    tntDoubleArraysSave(currentSavename, mpsEntanSpec1, mpsEntanSpec2, mpsEntanSpec3);
    tntDoubleArraysSave(currentSavename, mpsEntanEntr1, mpsEntanEntr2, mpsEntanEntr3);

    /* Calculate initial kinetic energy and enstrophy */
    tntMpsNsIncomp3DCalcEkinEnst(&(Ekin.vals[0]),&(Enst.vals[0]),u1,u2,u3,d1,d2,d3,chi,0);
    
    Ek0= tntMpsSelfProduct(u1,0) + tntMpsSelfProduct(u2,0) + tntMpsSelfProduct(u3,0);
    printf("\nTime evolution starting at t = %f. Ek0 = %.4f\n",t,Ek0); fflush(stdout);
    
    /* Set dynamic adjustments to disabled initially */
    sprintf(dynamicAdjustments,"None");
    
    /* Run RK2, i.e. midpoint method, to time-evolve */
    time_t start, end;
    double time_used;
    while(t<T*0.99999 || T == 0){
        start = time(NULL);
        tntWalltimeSet(0.0);
        
        /* First calculate advection-diffusion term at initial field and store in ad1 and ad2 arrays */
        tntMpsNsIncomp3DBuildBurgerTerms(&BurgMpoSum, &numBurgTerms,
                                     u1, u2, u3,
                                     d1, dd1,
                                     d2, dd2,
                                     d3, dd3,
                                     Re, dt/2, chi, chi_mpo, errp);
        
        /* First calculate advection-diffusion term at initial field and store in ad1 and ad2 arrays */
        tntMpsNsIncomp3D1sTimeStep(*v1, *v2, *v3,
                               u1, u2, u3, u1,u2,u3,
                               BurgMpoSum,
                               d1xd1, d1xd2, d1xd3, d2xd2, d2xd3, d3xd3,
                               &R1d1xd1R1_left, &R1d1xd1R1_right,
                               &R1d1xd2R2_left, &R1d1xd2R2_right,
                               &R1d1xd3R3_left, &R1d1xd3R3_right,
                               &R2d2xd2R2_left, &R2d2xd2R2_right,
                               &R2d2xd3R3_left, &R2d2xd3R3_right,
                               &R3d3xd3R3_left, &R3d3xd3R3_right,
                               numBurgTerms,
                               epsilon, dynamicAdjustments);
        
        /* Extract midpoint field */
        uMidpoint1 = tntNetworkCopy(*v1);
        uMidpoint2 = tntNetworkCopy(*v2);
        uMidpoint3 = tntNetworkCopy(*v3);
                
        /* Now delete old advection-diffusion terms and calculate new advection-diffusion term at the current field, i.e. at the midpoint */
        tntNetworkArrayFree(BurgMpoSum, numBurgTerms);
        tntMpsNsIncomp3DBuildBurgerTerms(&BurgMpoSum, &numBurgTerms,
                                         uMidpoint1, uMidpoint2, uMidpoint3,
                                         d1, dd1,
                                         d2, dd2,
                                         d3, dd3,
                                         Re, dt, chi, chi_mpo, errp);
        
        /* Use the updated advection-diffusion terms to time-evolve a full dt forward */
        tntMpsNsIncomp3D1sTimeStep(*v1, *v2, *v3,
                                   u1, u2, u3, uMidpoint1, uMidpoint2, uMidpoint3,
                                   BurgMpoSum,
                                   d1xd1, d1xd2, d1xd3, d2xd2, d2xd3, d3xd3,
                                   &R1d1xd1R1_left, &R1d1xd1R1_right,
                                   &R1d1xd2R2_left, &R1d1xd2R2_right,
                                   &R1d1xd3R3_left, &R1d1xd3R3_right,
                                   &R2d2xd2R2_left, &R2d2xd2R2_right,
                                   &R2d2xd3R3_left, &R2d2xd3R3_right,
                                   &R3d3xd3R3_left, &R3d3xd3R3_right,
                                   numBurgTerms,
                                   epsilon, dynamicAdjustments);
        end = time(NULL);
        time_used = ((double) (end - start));// / CLOCKS_PER_SEC;
        printf("\nTotal time spent doing one RK2 timestep: %.0fs. Of that, %.0fs was spent in tensor_reshape.",time_used, tntWalltimeGet()); fflush(stdout);
        /* Move one timestep forward */
        tntNetworkFree(&u1);
        tntNetworkFree(&u2);
        tntNetworkFree(&u3);
        tntNetworkFree(&uMidpoint1);
        tntNetworkFree(&uMidpoint2);
        tntNetworkFree(&uMidpoint3);
        u1 = tntNetworkCopy(*v1);
        u2 = tntNetworkCopy(*v2);
        u3 = tntNetworkCopy(*v3);
        
        /* timestep forward */
        t=t+dt;
        step++;
        
        /* Advection-diffusion terms now out of date, so must be deleted in preperation for next step */
        tntNetworkArrayFree(BurgMpoSum, numBurgTerms);
        
        /* Calculate kinetic energy and enstrophy at time t. */
        tntMpsNsIncomp3DCalcEkinEnst(&(Ekin.vals[step]),&(Enst.vals[step]),u1,u2,u3,d1,d2,d3,chi,0);
        
        Ek= tntMpsSelfProduct(u1,0) + tntMpsSelfProduct(u2,0) + tntMpsSelfProduct(u3,0);
        printf("\nRK2 step finished at t = %f. Ek1 = %.4f\n",t, Ek); fflush(stdout);
        /* If appropriate, sample the current MPS and return */
        
        if(step%stepsBetweenSampling == 0){
            
            /* Calculate current entanglement spectrum and entropy. */
            tntMpsCalcEntanSpec(u1, chi, 0, &mpsEntanSpec1, &mpsEntanEntr1);
            tntMpsCalcEntanSpec(u2, chi, 0, &mpsEntanSpec2, &mpsEntanEntr2);
            tntMpsCalcEntanSpec(u3, chi, 0, &mpsEntanSpec3, &mpsEntanEntr3);
            
            printf("Saving time-evolved fields at t=%f.\n",t);
            printf("\nNormalised Entanglement Entropies at every bond: ");
            for(int i = 0; i<L-1; i++){
                printf("[%.4f,%.4f,%.4f]_%i; ",mpsEntanEntr1.vals[i],mpsEntanEntr2.vals[i],mpsEntanEntr3.vals[i],i+1);
            }  printf("\n"); fflush(stdout);
            
            currentSavename[0] = '\0';
            strcat(currentSavename,savepath);
            sprintf(fieldname, "/NsIncomp3DFieldOutL%iC%iRe%.0fmt%.0f.mat",L,chi,Re,round(t*1000));
            strcat(currentSavename, fieldname );
            tntNetworksSave(currentSavename, u1, u2, u3);
            tntDoubleArraysSave(currentSavename, mpsEntanSpec1, mpsEntanSpec2, mpsEntanSpec3);
            tntDoubleArraysSave(currentSavename, mpsEntanEntr1, mpsEntanEntr2, mpsEntanEntr3);
            tntDoubleParamsSave(currentSavename, Re,t,T_sampling, lengthScale);
            
            tntDoubleArraysSave(currentSavename, Ekin, Enst);
            tntDoubleParamsSave(currentSavename, T,dt);
            
            /* Break & Exit if Ekin is zero */
            if(Ek/Ek0<1e-8 || Ek < 1e-6){
                printf("Kinetic energy is now practically zero, since Ek(t)/Ek(0) < 1e-10. Exiting.\n");
                break;
            }
            
        }
        
        /* Set dynamic adjustments to enabled now that the simulation has "warmed up" */
        sprintf(dynamicAdjustments,"Yes");
        
        /* If dt = 0, the time-evolution is finished. Get out of the loop. */
        if(t==0 || T==0){
            break;
        }
    }
    
    /* Save Kinetic energy and Enstrophy as functions of time */
    currentSavename[0] = '\0'; strcat(currentSavename,savepath);
    sprintf(txtfilename, "/NsIncomp3DEkinEnstL%iC%iRe%.0f.mat",L,chi,Re);
    strcat(currentSavename, txtfilename);
    printf("\nSaving final kinetic energy and enstrophy at t=%f to file %s.",t, txtfilename); fflush(stdout);
    tntDoubleArraysSave(currentSavename, Ekin, Enst);
    tntDoubleParamsSave(currentSavename, Re,T,dt, lengthScale);
    
    
    /* Free everything and exit */
    tntDoubleArrayFree(&Ekin);
    tntDoubleArrayFree(&Enst);
    tntDoubleArrayFree(&mpsEntanSpec1);
    tntDoubleArrayFree(&mpsEntanSpec2);
    tntDoubleArrayFree(&mpsEntanSpec3);
    tntDoubleArrayFree(&mpsEntanEntr1);
    tntDoubleArrayFree(&mpsEntanEntr2);
    tntDoubleArrayFree(&mpsEntanEntr3);
    
    tntNetworkFree(&u1);
    tntNetworkFree(&u2);
    tntNetworkFree(&u3);
    
    tntNetworkFree(&d1);
    tntNetworkFree(&d2);
    tntNetworkFree(&d3);
    
    tntNetworkFree(&dd1);
    tntNetworkFree(&dd2);
    tntNetworkFree(&dd3);
    
    tntNetworkFree(&d1xd1);
    tntNetworkFree(&d1xd2);
    tntNetworkFree(&d1xd3);
    tntNetworkFree(&d2xd2);
    tntNetworkFree(&d2xd3);
    tntNetworkFree(&d3xd3);
    
    tntNodeArrayFree(&R1d1xd1R1_left); tntNodeArrayFree(&R1d1xd1R1_right);
    tntNodeArrayFree(&R1d1xd2R2_left); tntNodeArrayFree(&R1d1xd2R2_right);
    tntNodeArrayFree(&R1d1xd3R3_left); tntNodeArrayFree(&R1d1xd3R3_right);
    tntNodeArrayFree(&R2d2xd2R2_left); tntNodeArrayFree(&R2d2xd2R2_right);
    tntNodeArrayFree(&R2d2xd3R3_left); tntNodeArrayFree(&R2d2xd3R3_right);
    tntNodeArrayFree(&R3d3xd3R3_left); tntNodeArrayFree(&R3d3xd3R3_right);
    
}
