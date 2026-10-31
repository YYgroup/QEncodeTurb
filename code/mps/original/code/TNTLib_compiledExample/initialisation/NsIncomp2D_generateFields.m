restoredefaultpath;
clear restoredefaultpath;

addpath('../../PDE_matfiles/MPSBuilderNik');
addpath('../../PDE_matfiles/MPSUtility');
addpath('../../PDE_matfiles/NS2DIncomp');
addpath('../../PDE_matfiles/NS2DIncomp/Field2DBuilder');

mkdir ../output/

%% Create field
L_orig = 6;%10;
L_new = L_orig;
chi0 = 2^L_orig;
kmax = 2; 
alpha = 1; %Must be an integer.
Umax = 1;%1/32;
BC1 = 'Periodic';
BC2 = 'Periodic';

epsilon = 10^-8;

%% Set timesteps
dt = 0.1*2^-(L_new-1);
T = 3;
T_sampling = 0.025;

%% Set penalty-factor for breach of incompressibility condition
dx = 2^-L_new;
lambda = dx^2*10^4;

[v1, v2, U, V] = buildFluidField2DMPS(Umax, L_orig, chi0, kmax,alpha);

%% Save fields
for chi = 33;%[13,28,60,149]
    for Re = 1000*kmax
        
        chi_mpo = -1;%max(chi, 32);
    end
    
    %
    % %% Visualise
    % %%VisualiseFlowField(U,V)
    %
%     %% Calculate Initial error
%     MPSU_chi = createMPS2D(U,chi);
%     MPSV_chi = createMPS2D(V,chi);
%     Ur = reconstructFromNik2D(MPSU_chi);
%     Vr = reconstructFromNik2D(MPSV_chi);
%     VisualiseFlowField(Ur,Vr,0)
%     %
%     Fidelity = 0.5*calcFidelity(Ur,U) + 0.5*calcFidelity(Vr,V);
%     fprintf("\nRanks of U and V are %.0f and %.0f.\nInfidelity = %.5e at Chi = %d \n\n",rank(U),rank(V),1-Fidelity,chi);
    %% Save
    filename = sprintf('NsIncomp2DField_L%iC%iRe%i',L_new,chi,Re);
    save(filename);
    
end