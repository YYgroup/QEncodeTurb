Author: Nikita Gourianov
Purpose: DEMO CODE that produces results similar to those found in "A Quantum Inspired Approach to Exploit Turbulence Structures (2021)" paper. 

This DEMO MATLAB code is used to analyse the results of the 2D and 3D DNS. It produces results analogous to those found in Fig. 1 of the main text and supplementary Figs. 1, 2 and 3 of the "A Quantum Inspired Approach to Exploit Turbulence Structures (2021)" paper. The only difference is that the system size is significantly smaller here. It should take a normal desktop computer ~10 minutes of computing time to run this demo.

Before the DEMO figures can be produced, DNS data must be generated. This done in the following two steps.

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

Note: Please use the code in the directory "AnalysingResultsFromDNS" to reproduce the data found in the article.