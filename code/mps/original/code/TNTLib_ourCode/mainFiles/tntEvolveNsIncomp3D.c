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
    
    tntNetwork v1, v2, v3;  /* The initial velocity fields along x and y in MPS form */
    tntNetwork v1_prol, v2_prol, v3_prol; /* post-prolongation to L_new */
    double t0, dt;  /* starting time and time-step */
    double T; /* Final time */
    double T_sampling; /* Sampling time at which to store output */
    double Re; /* Reynolds number */
    double lengthScale; /* Length of spatial domain */
    double lambda; /* Penalty coefficient */
    double epsilon; /* Sensitivity at which to decide the cost function minimisation has converged */
    double errp; /* dummy variable */
    int L_new, L, chi, chi_mpo; /* Requested length of MPS chain along with requested maximum bond dimension of wavefunction, and max bond dimension of mpos */
    double err; /* dummy variable */
    
    /* Initialize the TNT library */
    tntInitialize();
    
    /* Set random seed */
    tntSysRandSeedSet(0);
    
    /* Process the command line options */
    tntProcessCLOptions(argc, argv, NULL, initfile, savename, NULL, NULL);
    /* Load the requested chi */
    tntIntParamsLoad(initfile, L_new,  chi, chi_mpo);
    tntDoubleParamsLoad(initfile, Re, lambda, epsilon, lengthScale);
    tntDoubleParamsLoad(initfile, t0, dt, T, T_sampling);
    /* Load fields from initialisation file */
    tntNetworksLoad(initfile, v1);
    tntNetworksLoad(initfile, v2);
    tntNetworksLoad(initfile, v3);
    
    L = tntMpsLength(v1);
    
    /*************************************************************************************************/
    /*        Perform prolongation to get to L_new          */
    /*************************************************************************************************/
    
    if(L==L_new){
        v1_prol = v1;
        v2_prol = v2;
        v3_prol = v3;
    }
    else{
        v1_prol = tntNetworkCopy(v1);
        v2_prol = tntNetworkCopy(v2);
        v3_prol = tntNetworkCopy(v3);

        tntMpsPdeProlongate(&v1_prol, L_new, chi, "Periodic", &errp);
        tntMpsPdeProlongate(&v2_prol, L_new, chi, "Periodic", &errp);
        tntMpsPdeProlongate(&v3_prol, L_new, chi, "Periodic", &errp);
        tntNetworkFree(&v1);
        tntNetworkFree(&v2);
        tntNetworkFree(&v3);
    }
    
    /*************************************************************************************************/ 
    /*                      Time evolve                     */
    /*************************************************************************************************/
    
    tntMpsBoxTurbulence3DTimeEvolutionRK2(&v1_prol, &v2_prol, &v3_prol,t0, dt, T, T_sampling, savename, "Periodic", "Periodic", "Periodic", Re, lengthScale, lambda, epsilon, chi, chi_mpo, &errp);

    tntFinalize();
    
    return 0;
}
