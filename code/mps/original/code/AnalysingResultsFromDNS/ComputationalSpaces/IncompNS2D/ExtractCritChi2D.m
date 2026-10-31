% Finds critical Chis for representing the flow at a certain accuracy.

%% First add paths to access required functions and files
clear all; warning('off','all');
addpath('../../DNS/IncNs2DMatlab/');
addpath('../../DNS/IncNs2DMatlab/MPSBuilder/');

%% Set parameters
% %Option 1: Compare MPS representation against URDNS representation (some fourier-interpolated rough grid) at lengthscale l_crit
%l_crit = log2(362/2); %Note: N fourier modes is equivalent to a (2N)^2 grid, not a N^2 grid. 

%Option 2: Specify desired accuracy 
desiredAcc = [10^-2,10^-3,10^-4,10^-5,10^-6,10^-7,10^-8,10^-9,10^-10,10^-11,10^-12];

desiredAcc250 = desiredAcc;
desiredAcc750 = desiredAcc;
desiredAcc1250 = desiredAcc;
desiredAcc1750 = desiredAcc;


%% Load 2D at t=0.75, 1.25, 1.75
filenames{1} = sprintf('../../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/intermediate/Vortex2D_Re200000_L10.0_mdt0_mt250.mat');
filenames{2} = sprintf('../../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/intermediate/Vortex2D_Re200000_L10.0_mdt0_mt750.mat');
filenames{3} = sprintf('../../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/intermediate/Vortex2D_Re200000_L10.0_mdt0_mt1250.mat');
filenames{4} = sprintf('../../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/intermediate/Vortex2D_Re200000_L10.0_mdt0_mt1750.mat');

%t=0.25
filename = filenames{1};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    U = NaN*ones(N,N); V = NaN*ones(N,N);
else
    load(filename);
end
U250 = U; V250 =V;

%t=0.75
filename = filenames{2};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    U = NaN*ones(N,N); V = NaN*ones(N,N);
else
    load(filename);
end
U750 = U; V750 =V;

%t=1.25
filename = filenames{3};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    U = NaN*ones(N,N); V = NaN*ones(N,N);
else
    load(filename);
end
U1250 = U; V1250 =V;


%t=1.75
filename = filenames{4};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    U = NaN*ones(N,N); V = NaN*ones(N,N);
else
    load(filename);
end
U1750 = U; V1750 =V;

% %% Calculate the accuracy of URDNS representation 
% %t=0.25
% E250 = enspec2D_homo(U250,V250);
% desiredAcc250 = 1-(sum(E250(1:min(2^l_crit))));
% 
% %t=0.75
% E750 = enspec2D_homo(U750,V750);
% desiredAcc750 = 1-(sum(E750(1:min(2^l_crit))));
% 
% %t=1.25
% E1250 = enspec2D_homo(U1250,V1250);
% desiredAcc1250 = 1-(sum(E1250(1:2^l_crit)));
% 
% %t=1.75
% E1750 = enspec2D_homo(U1750,V1750);
% desiredAcc1750 = 1-(sum(E1750(1:min(2^l_crit))));
% 
%% Find critical chi for each of these
%t=0.25
[chi250,entanEntrU250,entanEntrV250,SingValsU250,SingValsV250] = findCritChi2D(U250,V250,desiredAcc250);

%t=0.75
[chi750,entanEntrU750,entanEntrV750,SingValsU750,SingValsV750] = findCritChi2D(U750,V750,desiredAcc750);

%t=1.25
[chi1250,entanEntrU1250,entanEntrV1250,SingValsU1250,SingValsV1250] = findCritChi2D(U1250,V1250,desiredAcc1250);

%t=1.75
[chi1750,entanEntrU1750,entanEntrV1750,SingValsU1750,SingValsV1750] = findCritChi2D(U1750,V1750,desiredAcc1750);

%% Calculate chiMax
N = length(U(1,:));
L = log2(N);
if(L ~= round(L))
    fprintf('\nError: (U,V) must be a 2D field where both dimensions are powers of 2. Exiting.\n');
    exit(1);
end
chiMax = zeros(L-1,1);

%Calculate it
for l = 1:(L-1)
    chiMax(l) = min(4^l,4^(L-l));
end

clear KX KY U U250 U750 U1250 U1750 V V250 V750 V1250 V1750

save('MWS_e2to12_out');