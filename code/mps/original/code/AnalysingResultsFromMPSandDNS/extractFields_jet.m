%% Extracts flow fields for visualisation

%% First add paths to access required functions and files
clear all; warning('off','all');

addpath('../PDE_matfiles/MPSUtility');
addpath('../PDE_matfiles/MPSBuilderNik');
addpath('../PDE_matfiles/NS2DIncomp');
addpath('../PDE_matfiles/General');

TNTpath_incompNS = '../../data/IncompNs2D/';
MATpath_incompNS = '../AnalysingResultsFromDNS/DNS/IncNs2DMatlab/2DFlowAcrossRe/output/intermediate/';
MATfstart = 'Vortex2D_Re200000_';
TNTInsFstart_ex = 'NsIncomp2DFieldOut';

%% First extract DNS data at Re = 200,000 (equivalent to Re=1000 in manuscript) which has been discretised on 1024^2, 362^2, 256^2 and 128^2 grids at times t = 0.25, 0.75, 1.25 and 1.75.

%128^2
%t=0.25
MATfend = sprintf('L7.0_mdt2_mt250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U250L7 =U;
V250L7 =V;

% t=0.75
MATfend = sprintf('L7.0_mdt2_mt750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U750L7 =U;
V750L7 =V;

%t=1.25
MATfend = sprintf('L7.0_mdt2_mt1250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1250L7 =U;
V1250L7 =V;

%t=1.75
MATfend = sprintf('L7.0_mdt2_mt1750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1750L7 =U;
V1750L7 =V;

% 256^2
% t=0.25
MATfend = sprintf('L8.0_mdt1_mt250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U250L8 =U;
V250L8 =V;

% t=0.75
MATfend = sprintf('L8.0_mdt1_mt750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U750L8 =U;
V750L8 =V;

% t=1.25
MATfend = sprintf('L8.0_mdt1_mt1250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1250L8 =U;
V1250L8 =V;

% t=1.75
MATfend = sprintf('L8.0_mdt1_mt1750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1750L8 =U;
V1750L8 =V;

% 362^2
% t=0.25
MATfend = sprintf('L8.5_mdt1_mt250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U250L8p5 =U;
V250L8p5 =V;

% t=0.75
MATfend = sprintf('L8.5_mdt1_mt750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U750L8p5 =U;
V750L8p5 =V;

% t=1.25
MATfend = sprintf('L8.5_mdt1_mt1250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1250L8p5 =U;
V1250L8p5 =V;

% t=1.75
MATfend = sprintf('L8.5_mdt1_mt1750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1750L8p5 =U;
V1750L8p5 =V;

% 1024^2
% t=0.25
MATfend = sprintf('L10.0_mdt0_mt250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U250 =U;
V250 =V;

% t=0.75
MATfend = sprintf('L10.0_mdt0_mt750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U750 =U;
V750 =V;

% t=1.25
MATfend = sprintf('L10.0_mdt0_mt1250.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1250 =U;
V1250 =V;

% t=1.75
MATfend = sprintf('L10.0_mdt0_mt1750.mat');
filename = strcat(MATpath_incompNS, MATfstart, MATfend);
load(filename);
U1750 =U;
V1750 =V;

%% Now extract TNT data at Re = 200,000 (equivalent to Re=1000 in manuscript) discretised on 1024^2 grids, for bond dimensions chi = 33, 74, 118, at times t = 0.25, 0.75, 1.25 and 1.75.

%chi=33
%t=0.25
TNTfend = sprintf('L10C33Re200000mt250.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc33t250, Vc33t250 ]= reconstructFields(u1,u2);

%t=0.75
TNTfend = sprintf('L10C33Re200000mt750.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc33t750, Vc33t750 ]= reconstructFields(u1,u2);

%t=1.25
TNTfend = sprintf('L10C33Re200000mt1250.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc33t1250, Vc33t1250 ]= reconstructFields(u1,u2);

%t=1.75
TNTfend = sprintf('L10C33Re200000mt1750.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc33t1750, Vc33t1750 ]= reconstructFields(u1,u2);

%chi=74
%t=0.25
TNTfend = sprintf('L10C74Re200000mt250.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc74t250, Vc74t250 ]= reconstructFields(u1,u2);

%t=0.75
TNTfend = sprintf('L10C74Re200000mt750.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc74t750, Vc74t750 ]= reconstructFields(u1,u2);

%t=1.25
TNTfend = sprintf('L10C74Re200000mt1250.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc74t1250, Vc74t1250 ]= reconstructFields(u1,u2);

%t=1.75
TNTfend = sprintf('L10C74Re200000mt1750.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc74t1750, Vc74t1750 ]= reconstructFields(u1,u2);

%chi=118
%t=0.25
TNTfend = sprintf('L10C118Re200000mt250.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc118t250, Vc118t250 ]= reconstructFields(u1,u2);

%t=0.75
TNTfend = sprintf('L10C118Re200000mt750.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc118t750, Vc118t750 ]= reconstructFields(u1,u2);

%t=1.25
TNTfend = sprintf('L10C118Re200000mt1250.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc118t1250, Vc118t1250 ]= reconstructFields(u1,u2);

%t=1.75
TNTfend = sprintf('L10C118Re200000mt1750.mat');
filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
load(filename);
[Uc118t1750, Vc118t1750 ]= reconstructFields(u1,u2);

save('extractedFields.mat');
