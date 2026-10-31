Author: Nikita Gourianov
Purpose: Demo code for the "A Quantum Inspired Approach to Exploit Turbulence Structures (2021)" paper. 

The "Tensor Network Theory Library" (TNT library) was used to perform the MPS simulations found in Figs. 2 and 3. The TNT library is however configured to operate on HPC systems and is not open source. We therefore provide a compiled version of our MPS code in this directory. It is configured to run on a "Mid 2015 Retina Macbook Pro (2,5 GHz Quad-Core Intel Core i7) using macOS 11.4 and Xcode versioned at 12.5.1 (build version: 12E507)". Before running it, make sure you are in a computational environment analogous to this Macbook Pro, especially in terms of both the processor architecture and Kernel. 

To run this compiled code, do the following steps:

Step 0. Ensure your computing architecture is compatible with the above Macbook Pro.
Step 1. Navigate to this directory in the terminal, and run "source setEnviron.sh".
Step 2. Run the initialisation/NsIncomp2D_generateFields.m MATLAB file to initialise the data necessary for a Re=10 simulation on a 64^2 grid with chi = 33.
Step 3. Navigate to the /bin/ folder, and run ./tntEvolveNsIncomp2D -o ../initialisation/NsIncomp2DField_L6C33Re2000.mat -d ../output/. This will begin the simulation. Time-evolved flow fields will be outputed to the ../output/ directory at the sampling time period defined in the initialisation/NsIncomp2D_generateFields.m file.
Step 4. Navigate to the /output/ directory in Matlab and run "addpath('../../PDE_matfiles/MPSBuilderNik');addpath('../../PDE_matfiles/MPSUtility');addpath('../../PDE_matfiles/NS2DIncomp');addpath('../../PDE_matfiles/NS2DIncomp/Field2DBuilder');". Then, open one of the .mat files in the output (mtX denotes the time such that t = X/1000) and decompress the MPS using "[U,V] = reconstructFields(u1,u2)". Inspect U and V, which are the streamwise and spanwise velocity components, respectively. 