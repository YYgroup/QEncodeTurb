% Finds critical Chis for representing the flow at a certain accuracy.

%% First add paths to access required functions and files
clear all; warning('off','all');
addpath('../../DNS/IncNs3DMatlab/');
addpath('../../DNS/IncNs3DMatlab/MPSBuilder/');

%% Set parameters
%Option 1: Compare MPS representation against URDNS representation (some fourier-interpolated rough grid) at lengthscale l_crit
% N = 16/2; %Note: N fourier modes is equivalent to a (2N)^2 grid, not a N^2 grid. 

% %Option 2: Specify desired accuracy 
desiredAcc = [10^-2,10^-3,10^-4,10^-5,10^-6,10^-7,10^-8,10^-9,10^-10,10^-11,10^-12];

desiredAcc200 = desiredAcc;
desiredAcc800 = desiredAcc;
desiredAcc1400 = desiredAcc;
desiredAcc2000 = desiredAcc;

%% Load 2D at t=0.2, 0.8, 1.4, 2
filenames{1} = sprintf('../../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/intermediate/Vortex3D_Re5027_L8.0_mdt2_mt200.mat');
filenames{2} = sprintf('../../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/intermediate/Vortex3D_Re5027_L8.0_mdt2_mt800.mat');
filenames{3} = sprintf('../../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/intermediate/Vortex3D_Re5027_L8.0_mdt2_mt1400.mat');
filenames{4} = sprintf('../../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/final/Vortex3D_Re5027_L8.0_mdt2_mt2000.mat');

%t=0.2
filename = filenames{1};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    exit(1);
else
    load(filename);
end
U200 = U; V200 =V; W200 =W;

%t=0.8
filename = filenames{2};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    exit(1);
else
    load(filename);
end
U800 = U; V800 =V; W800 =W;

%t=1.4
filename = filenames{3};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    exit(1);
else
    load(filename);
end
U1400 = U; V1400 =V; W1400 =W;


%t=2
filename = filenames{4};
if(fopen(filename) == -1)
    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
    exit(1);
else
    load(filename);
end
U2000 = U; V2000 =V; W2000 =W;

%% Calculate the accuracy of URDNS representation 
% %t=0.2
% [~,E200] = enspec3D_homo(U200,V200,W200);
% desiredAcc200 = 1-(sum(E200(1:N)));
% 
% %t=0.8
% [~,E800] = enspec3D_homo(U800,V800,W800);
% desiredAcc800 = 1-(sum(E800(1:N)));
% 
% %t=1.4
% [~,E1400] = enspec3D_homo(U1400,V1400,W1400);
% desiredAcc1400 = 1-(sum(E1400(1:N)));
% 
% %t=2
% [~,E2000] = enspec3D_homo(U2000,V2000,W2000);
% desiredAcc2000 = 1-(sum(E2000(1:N)));

%% Find critical chi for each of these
%
N = length(U(1,:,:));
L = log2(N);
if(L ~= round(L))
    fprintf('\nError: (U,V,W) must be a 3D field where both dimensions are powers of 2. Exiting.\n');
    exit(1);
end

if( sum(desiredAcc200 < 10^-15) || sum(desiredAcc800 <10^-15) || sum(desiredAcc1400 <10^-15) || sum(desiredAcc2000 <10^-15) )
    fprintf('\nError: Requested accuracy is too high. Try again. Exiting. \n');
    exit(1);
end

%t=0.2
[chi200,entanEntrU200,entanEntrV200,entanEntrW200,SingValsU200,SingValsV200,SingValsW200] = findCritChi3D(U200,V200,W200,desiredAcc200);

%t=0.8
[chi800,entanEntrU800,entanEntrV800,entanEntrW800,SingValsU800,SingValsV800,SingValsW800] = findCritChi3D(U800,V800,W800,desiredAcc800);

%t=1.4
[chi1400,entanEntrU1400,entanEntrV1400,entanEntrW1400,SingValsU1400,SingValsV1400,SingValsW1400] = findCritChi3D(U1400,V1400,W1400,desiredAcc1400);

%t=2
[chi2000,entanEntrU2000,entanEntrV2000,entanEntrW2000,SingValsU2000,SingValsV2000,SingValsW2000] = findCritChi3D(U2000,V2000,W2000,desiredAcc2000);

%% Calculate chiMax
chiMax = zeros(L-1,1);

%Calculate it
for l = 1:(L-1)
    chiMax(l) = min(8^l,8^(L-l));
end

%% Clear and save
clear KX KY KZ U V W U200 U800 U1400 U2000 V200 V800 V1400 V2000 W200 W800 W1400 W2000;


save('MWS_e2to12_out');