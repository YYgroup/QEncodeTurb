% This is the demo script which reproduces Figure 1 for a smaller system than that found in the "A Quantum Inspired Apprach To Exploit Turbulence Structures" paper. It follows the instructions of the README file from "AnalysingResultsFromDNS_demo".


cd AnalysingResultsFromDNS_demo/DNS/IncNs2DMatlab/2DFlowAcrossRe/; IncNsReVarious;
cd ../../IncNs3DMatlab/3DFlowAcrossRe/; IncNsReVarious;
cd ../../../ComputationalSpaces/; get_d99;
cd IncompNS2D/; ExtractCritChi2D;
cd ../IncompNS3D/; ExtractCritChi3D;
cd ../; Plot_CompSpaces; EntanSpec_2D; EntanSpec_3D; Plot_EntanEntr;