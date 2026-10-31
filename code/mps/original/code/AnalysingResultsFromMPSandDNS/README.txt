Author: Nikita Gourianov
Purpose: Code needed to reproduce some of the results from the "A Quantum Inspired Approach to Exploit Turbulence Structures (2021)" paper. 

The MATLAB code required to reproduce Figs. 2 and 3. is provided here. This MATLAB code inputs the data from the ../../data/ and ../AnalysingResultsFromDNS/DNS/IncNs2DMatlab/2DFlowAcrossRe/output/ and ../AnalysingResultsFromDNS/DNS/IncNs3DMatlab/3DFlowAcrossRe/output/ folders, processes it, and produces Figs. 2 and 3b. The visual objects of Fig. 3a are also produced, but the software "paraview" must be used to inspect them. 

To reproduce Fig. 2., do the following.

Step 1. Go to ../AnalysingResultsFromDNS/DNS/IncNs2DMatlab/2DFlowAcrossRe/ and run IncNsRe1000.m.
Step 2. Come back to this directory, and run extractFields_jet.m
Step 3. Then run processFields_jet.m
Step 4. Run MakeGrandPlot_2D.m. This will reproduce Fig. 2, comparing the MPS results against DNS and under-resolved DNS. The figure is stored in ../results/.

Let us move on now and reproduce Fig.3b and the paraview objects of Fig. 3a.

Step 5. Go to ../AnalysingResultsFromDNS/DNS/IncNs3DMatlab/3DFlowAcrossRe/ and run IncNsRe800.m.
Step 6. Run extractData_3Dhomo.m. This will take a few hours.
Step 7. Run MakeGrandPlot_3D.m. This will generate Fig. 3b and store it in the ../results/ directory.
Step 8. Finally, run ProduceVTKs.m. This will generate the .vtk files of Fig. 3a and put them into the ../results/VTKfiles/ directory. These files can be opened and inspected with the open source software "Paraview". 
