% This is an example code that loads simulation data from a '.mat' file, 
% reconstructs velocity fields, and visualizes the flow field 
% using a the function reconstructFields.

filename = 'NsIncomp2DFieldOutL10C118Re200000mt1750.mat';
dir = fullfile('data', 'IncompNs2D', filename);
load(dir);

[data_U, data_V] = reconstructFields(u1,u2);
VisualiseFlowField_special(data_U,data_V,0,[1, 1, 0],12,1,1,1)