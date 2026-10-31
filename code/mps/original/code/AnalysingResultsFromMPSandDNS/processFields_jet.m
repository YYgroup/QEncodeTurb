%% Compares simulation results of TNT and Matlab 

%% First add paths to access required functions and files
clear all; warning('off','all');

addpath('../PDE_matfiles/MPSUtility');
addpath('../PDE_matfiles/MPSBuilderNik');
addpath('../PDE_matfiles/NS2DIncomp');
addpath('../PDE_matfiles/General');

TNTpath_incompNS = '../../data/IncompNs2D/';
MATpath_incompNS = '../AnalysingResultsFromDNS/DNS/IncNs2DMatlab/2DFlowAcrossRe/output/';

%% Define what to extract
LtoAn = [7 8 8.5 10];
Ts = 0.025;
Tf = 2.0;
mdt = [2 1 1 0];%1000*0.1*2^-9;
ReToAn = 200*1000;
chiToAn = [33, 74, 118];

MATfstart = 'Vortex2D_Re';
TNTInsFstart_ex = 'NsIncomp2DFieldOut';
TNTInsFstart_imex = 'NsIncompImex2DFieldOut';

%% Initialise necessary vectors
%Incompressible Navier-Stokes
%Data-format: [L, chi, data]
%chi \in chiToAn. L \in LtoAn.
N = round(2^max(LtoAn));

% insFid_TNT = NaN*ones(length(LtoAn),length(ReToAn),length(chiToAn),round(Tf/Ts)+1);
% insK_TNT = NaN*ones(N,length(LtoAn),length(ReToAn),length(chiToAn),round(Tf/Ts)+1);
% insK_MAT = NaN*ones(N,length(LtoAn),length(ReToAn),round(Tf/Ts)+1);
% insP1_TNT = insK_TNT; insP1_MAT = insK_MAT;
% insP2_TNT = insK_TNT; insP2_MAT = insK_MAT;
% insO_TNT = insK_TNT; insO_MAT = insK_MAT;
% insW_TNT = insK_TNT; insW_MAT = insK_MAT;
% insE_TNT = insK_TNT; insE_MAT = insK_MAT;
% insDiv_TNT = NaN*ones(length(LtoAn),length(ReToAn),length(chiToAn),round(Tf/Ts)+1);
% insDiv_MAT = NaN*ones(length(LtoAn),length(ReToAn),round(Tf/Ts)+1);
% insVT_TNT = NaN*ones(length(LtoAn),length(ReToAn),length(chiToAn),round(Tf/Ts)+1);
% insVT_MAT = NaN*ones(length(LtoAn),length(ReToAn),round(Tf/Ts)+1);
% insdP1dy_TNT = insK_TNT; insdP1dy_MAT = insK_MAT;
% insReStress11_TNT = insK_TNT; insReStress11_MAT = insK_MAT;
 insReStress12_TNT =  NaN*ones(N,length(LtoAn),length(ReToAn),length(chiToAn),round(Tf/Ts)+1); 
 insReStress12_MAT = NaN*ones(N,length(LtoAn),length(ReToAn),round(Tf/Ts)+1);
% insReStress22_TNT = insK_TNT; insReStress22_MAT = insK_MAT;

%insScalarFlux_MAT = insK_MAT;
%insST_MAT = insVT_MAT;

%% Loop through and extract error along with the aggregate variables
fprintf('\n========Starting loading the Incompressible Navier-Stokes data.========\n');

l = 0;
j = 0;
for L = LtoAn
    l = l + 1;
    N = 2^LtoAn(l);
    
    ReNum = 0;
    
    for Re = ReToAn
        ReNum = ReNum + 1;
        
        c = 0;
        for chi = chiToAn
            c = c + 1;
            
            mt = 0;
            i = 0;
            while(mt<=Tf*1000)
                i=i+1;
                
                %display current state
                fprintf('\nCurrently at L=%s,Re=%d,Chi=%d,mt=%d',num2str(L),Re,chi,round(mt));
                
                %% Calculate MATLAB quantities
                if(mt<Tf*1000)
                   MATpath_incompNS2= strcat(MATpath_incompNS,'intermediate/');
                else
                   MATpath_incompNS2= strcat(MATpath_incompNS,'final/');
                end
                MATfend = sprintf('%d_L%.1f_mdt%d_mt%d.mat',Re,(L),round(mdt(l)),round(mt));
                filename = strcat(MATpath_incompNS2, MATfstart, MATfend);
                
                if(fopen(filename) == -1)
                    fprintf('. Breaking as Matlab file with filename "%s" doesn''t exist.',filename);
                    break;
                    U = NaN*ones(N,N); V = NaN*ones(N,N);
                else
                    load(filename);
                end
                
                % Calculate quantities                
                [~, ~, ~, ~, ~, ~, ~, ~,~,ReStress12,~] = calcAggAlongColsQuantities2D(U,V);
                [reStress, scalarThick, scalarFlux] = calcAggAlongColsConsScQuantities2D(Phi, U, V);
                
                % Store them in vector
%                 insP1_MAT(1:N,l,ReNum,i) = P(1);
%                 insP2_MAT(1:N,l,ReNum,i) = P(2);
%                 insK_MAT(1:N,l,ReNum,i) = K;
%                 insO_MAT(1:N,l,ReNum,i) = O;
%                 insW_MAT(1:N,l,ReNum,i) = W;
%                 insE_MAT(1:N,l,ReNum,i) = E;
%                 insDiv_MAT(l,ReNum,i) = div;
%                 insVT_MAT(l,ReNum,i) = VortThick;
%                 insdP1dy_MAT(1:N,l,ReNum,i) = dP1dy;
%                 insReStress11_MAT(1:N,l,ReNum,i) = ReStress11;
                insReStress12_MAT(1:N,l,ReNum,i) = ReStress12;
%                insReStress22_MAT(1:N,l,ReNum,i) = ReStress22;
                
                %
              %  insScalarFlux_MAT(1:N,l,ReNum,i) = scalarFlux;
              %  insST_MAT(l,ReNum,i) = scalarThick;
                
                % Calculate TNT quantities
                TNTfend = sprintf('L%sC%dRe%dmt%d.mat',num2str(L),chi,Re,round(mt));
                filename = strcat(TNTpath_incompNS, TNTInsFstart_ex, TNTfend);
                
                if(fopen(filename) == -1)
                    filename = strcat(TNTpath_incompNS, TNTInsFstart_imex, TNTfend);
                    if(fopen(filename) == -1)
                    %    fprintf('.\nBreaking as no TNT file with filename \n"%s" nor "%s"\n exists.',filename,strcat(TNTpath_incompNS, TNTInsFstart_imex, TNTfend));
                    %break;
                    end
                else
                load(filename);
                
                % Calculate quantities
                [Ur, Vr] = reconstructFields(u1,u2);
               % fidelity = 0.5*calcFidelity(Ur,U) + 0.5*calcFidelity(Vr,V);
                [~, ~, ~, ~, ~, ~, ~, ~,~, ReStress12,~] = calcAggAlongColsQuantities2D(Ur,Vr);
                
                % Store them in vector
              %  insFid_TNT(l,ReNum,c,i) = fidelity;                
              %  insP1_TNT(:,l,ReNum,c, i) = P(1);
              %  insP2_TNT(:,l,ReNum,c, i) = P(2);
              %  insK_TNT(:,l,ReNum,c, i) = K;
              %  insO_TNT(:,l,ReNum,c, i) = O;
              %  insW_TNT(:,l,ReNum,c, i) = W;
              %  insE_TNT(:,l,ReNum,c, i) = E;
              %  insDiv_TNT(l,ReNum,c, i) = div;
              %  insVT_TNT(l,ReNum,c, i) = VortThick;
              %  insdP1dy_TNT(:,l,ReNum,c, i) = dP1dy;
              %  insReStress11_TNT(:,l,ReNum,c, i) = ReStress11;
                insReStress12_TNT(:,l,ReNum,c, i) = ReStress12;
              %  insReStress22_TNT(:,l,ReNum,c, i) = ReStress22;
                end
                %% Continue timestepping
                
                % If at last TNT file, get the "Enstrophy vs Ekin" array
                mt = mt+1000*Ts;
            end
            
        end
    end   
end

fclose('all');
clear U; clear V; clear Ur; clear Vr; clear ekin; clear enst;
fprintf('\n================Done. Saving and exiting.================\n');

save('processedFields');
