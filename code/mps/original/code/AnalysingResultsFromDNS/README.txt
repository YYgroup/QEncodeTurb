Author: Nikita Gourianov
Purpose: Code needed to reproduce some of the results within the "A Quantum Inspired Approach to Exploit Turbulence Structures (2021)" paper. 


The MATLAB code used to analyse the results of the 2D and 3D DNS is provided here. Everything that's necessary to reproduce Fig. 1 of the main text and supplementary Figs. 1, 2 and 3 is contained in this directory. Furthermore, the function memComp.m is also provided. It calculates the compression ratios referenced in Figs. 2 and 3.

Before the figures can be reproduced though, the DNS data must be generated. This done in the following two steps.

Step 1a. Generate data for the 2D flow. Go to /DNS/IncNs2DMatlab/2DFlowAcrossRe/ and run IncNsReVarious.m. 
Step 1b. Generate data for the 3D flow. Go to /DNS/IncNs3DMatlab/3DFlowAcrossRe/ and run IncNsReVarious.m.

Now let us reproduce Fig. 1 of the main text and the supplementary figures.

Step 2. Go to /ComputationalSpaces/ and run get_d99.m.
Step 3a. Go to /ComputationalSpaces/IncompNS2D and run ExtractCritChi2D.m.
Step 3b. Go to /ComputationalSpaces/IncompNS3D and run ExtractCritChi3D.m.
Step 4. Go to /ComputationalSpaces/ and run Plot_CompSpaces.m. This will produce Figure 1 (except Fig. 1a, but it is merely an illustration and thus not relevant here). 
Step 5. Go to /ComputationalSpaces/ and run EntanSpec_2D.m. This will produce Supplementary Figure 1.
Step 6. Go to /ComputationalSpaces/ and run EntanSpec_3D.m. This will produce Supplementary Figure 2.
Step 7. Go to /ComputationalSpaces/ and run Plot_EntanEntr.m This will produce Supplementary Figure 3.

All figures are stored in ../results.

Finally, memComp.m exists in the main directory. Open it, read the instructions, and run it to recover the compression ratios referenced around Figs. 2 and 3.

Note: these computations are somewhat demanding and will take about two day on a good desktop computer. If a lighter "demo" example is preferred, please see the "demo" version of this software in the directory "AnalysingResultsFromDNS_demo". 