Author: Nikita Gourianov
Purpose: Code needed to reproduce the results within the "A Quantum Inspired Approach to Exploit Turbulence Structures (2021)" paper. 

The MATLAB code used to analyse the results of the 2D and 3D DNS is provided in the /AnalysingResultsFromDNS/ directory. Everything that's necessary to reproduce Fig. 1 of the main text and supplementary Figs. 1, 2 and 3 is contained there, along with the function memComp.m which calculates the compression ratios referenced in Figs. 2 and 3. AnalysingResultsFromDNS_demo contains the same files, but configured to run a smaller "demo" example. This demo example should only take ~5 minutes to run on a standard desktop machine. It can be run by either following the instructions in the README in /AnalysingResultsFromDNS_demo/, or by running the "runDemo.m" script within this folder.

The "Tensor Network Theory Library" (TNT library) was used to perform the MPS simulations of Figs. 2 and 3. The TNT library is configured to operate on HPC systems and is not open source. The code we developed for this paper is provided in the the /TNTLib_ourCode/ directory and employs the TNT library version 2.0.1. We provide a compiled demo-version of our tensor network code in the /TNTLib_compiledExample/ directory. This compiled code contains the TNT library and all other dependencies. This compiled bin file is configured to run on a "Mid 2015 Retina Macbook Pro using macOS 11.4 and Xcode versioned at 12.5.1 (build version: 12E507).

The raw data from our tensor network simulations are stored in the ../data/ folder and the MATLAB code required to reproduce Figs. 2 and 3. from this data is provided in the /AnalysingResultsFromMPSandDNS/ folder.

All figures produced by this code is put into the /results/ folder.

NOTE: The variables in this code differ from how they are called in the paper.
* In 2D flows, the Reynolds number as defined in the code is equivalent to the 2D Reynolds number of the manuscript multiplied by 200.
* In 3D flows, the Reynolds number as defined in the code is equivalent to the 3D Reynolds number of the manuscript multiplied by 2*pi.
* In 2D flows, 1/kmax (=1/200) as defined in the code is equivalent to h (=1/200) in the manuscript.