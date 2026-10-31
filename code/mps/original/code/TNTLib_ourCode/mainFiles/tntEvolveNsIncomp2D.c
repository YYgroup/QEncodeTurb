/*
Authors: Nikita Gourianov and Dr. Michael Lubasch
Date:    $LastChangedDate$
(c) Nikita Gourianov 2019
*/

#include "tntMps.h"
#include "tntPde.h"
#include "tntIncompNs.h"

int main(int argc, char **argv)
{

    char initfile[TNT_STRLEN];          /* Path to initialisation data */
    char savename[TNT_STRLEN];          /* Path to output file */
    
    tntNetwork v1, v2;  /* The initial velocity fields along x and y in MPS form */
    tntNetwork v1_prol, v2_prol; /* post-prolongation to L_new */
    double dt; /* time-step */
    double T; /* Final time */
    double T_sampling; /* Sampling time at which to store output */
    double Re;  /* Reynolds number */
    double lambda; /* Penalty coefficient */
    double epsilon; /* Sensitivity at which to decide the cost function minimisation has converged */
    double errp; /* dummy variable */
    int L_new, L, chi, chi_mpo; /* Requested length of MPS chain along with requested maximum bond dimension */
    double err; /* dummy variable */
    char BC1[50], BC2[50]; /* boundary conditions of x and y walls */
    
    /* Initialize the TNT library */
    tntInitialize();
    
    /* Set random seed */
    tntSysRandSeedSet(0);
    
    /* Process the command line options */
    tntProcessCLOptions(argc, argv, NULL, initfile, savename, NULL, NULL);
    /* Load the requested chi */
    tntIntParamsLoad(initfile, L_new,  chi, chi_mpo);
    tntStringsLoad(initfile,  BC1, BC2);
    tntDoubleParamsLoad(initfile, Re, lambda, epsilon);
    tntDoubleParamsLoad(initfile, dt, T, T_sampling);
    /* Load fields from initialisation file */
    tntNetworksLoad(initfile, v1);
    tntNetworksLoad(initfile, v2);
    
    L = tntMpsLength(v1);
    
    /*************************************************************************************************/
    /*        Perform prolongation to get to L_new          */
    /*************************************************************************************************/
    
    if(L==L_new){
        v1_prol = v1;
        v2_prol = v2;
    }
    else{
        v1_prol = tntNetworkCopy(v1);
        v2_prol = tntNetworkCopy(v2);
        tntMpsPdeProlongate(&v1_prol, L_new, chi, BC1, &errp);
        tntMpsPdeProlongate(&v2_prol, L_new, chi, BC2, &errp);
        tntNetworkFree(&v1);
        tntNetworkFree(&v2);
    }
    
    /*************************************************************************************************/ 
    /*                      Time evolve                     */
    /*************************************************************************************************/ 
    tntMpsBoxTurbulence2DTimeEvolutionRK2(&v1_prol, &v2_prol, dt, T, T_sampling, savename, BC1, BC2, Re, lambda, epsilon, chi, chi_mpo, &errp);
    tntFinalize();
    
    return 0;
}
