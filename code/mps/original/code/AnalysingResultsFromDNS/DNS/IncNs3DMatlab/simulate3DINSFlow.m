function [U_SingMatMax,V_SingMatMax, W_SingMatMax] = simulate3DINSFlow(lengthScale,kmax, Re,dt,t_final,TbetweenDDs,outputPath)
%Simulates 3D flow using the scripts in the 3DTurbulentFlowProducer folder
%while simultaneously converting flow fields into matrix product style 
%tensor networks (TNs) and extracting and analysing the entanglement spectra.
%t_final, dt and numMakeTN must be integer powers of 10. numMakeTN must be
%less than simutimeStep.

%Define time-evolution variables
N               = round(2^lengthScale); % N^3 grid, table allows looping over a mesh sensitivity test
simutimeSteps   = ceil(t_final/dt);
StepsBetweenDDs = floor(TbetweenDDs/dt);
desiredAcc = [10^-4];

% Non-dimensionalise the problem
Umax=1;
L=1;

%Initialise field
CreateFields3D;

%%
%Initialise time-series arrays
Ekin=zeros(simutimeSteps+1,1); 
vortSq=Ekin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes simulation begins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if(lengthScale == round(lengthScale))
    %
    U_SingMatMax = zeros(2,floor(simutimeSteps/StepsBetweenDDs));
    V_SingMatMax = zeros(2,floor(simutimeSteps/StepsBetweenDDs));
    W_SingMatMax = zeros(2,floor(simutimeSteps/StepsBetweenDDs));
    %
    D_max = 8^(floor(lengthScale/2)); %Max sing val dimension is sqrt(8^lengthScale) if lengthScale even, and sqrt(8^(lengthScale-1)) otherwise.
    k=1;
    %
    U_errT = zeros(D_max,floor(simutimeSteps/StepsBetweenDDs));
    U_errSmax = U_errT;
    V_errT = U_errT;
    V_errSmax = U_errT;
    W_errT = U_errT;
    W_errSmax = U_errT;
    %
    desiredAcc = [10^-4];
    chis = zeros(1,floor(simutimeSteps/StepsBetweenDDs));
end

    k=1;
%
for step = 0:simutimeSteps
    t = step*dt;
    
    %% 
    %Data dumps take a long time, so do them sparingly. 
    if(mod(step,StepsBetweenDDs)==0 || step == simutimeSteps)
        %Represent U in TN form and analyse if lengthScale permits
        if(lengthScale == round(lengthScale))
            %Represent U in TN form and analyse
            [SingMats,eEs] = createMPS3D(U,D_max);
            [Smax,bond,errSmax_vec,errT_vec] = TNanalyseFlow(SingMats,eEs);
            %
            U_SingMatMax(1,k) = Smax;
            U_SingMatMax(2,k) = bond;
            %
            U_errT(:,k) = errT_vec;
            U_errSmax(:,k) = errSmax_vec;
            
            %Represent V in TN form and analyse
            [SingMats,eEs] = createMPS3D(V,D_max);
            [Smax,bond,errSmax_vec,errT_vec] = TNanalyseFlow(SingMats,eEs);
            %
            V_SingMatMax(1,k) = Smax;
            V_SingMatMax(2,k) = bond;
            %
            V_errT(:,k) = errT_vec;
            V_errSmax(:,k) = errSmax_vec;
            
            %Represent W in TN form and analyse
            [SingMats,eEs] = createMPS3D(W,D_max);
            [Smax,bond,errSmax_vec,errT_vec] = TNanalyseFlow(SingMats,eEs);
            %
            W_SingMatMax(1,k) = Smax;
            W_SingMatMax(2,k) = bond;
            %
            W_errT(:,k) = errT_vec;
            W_errSmax(:,k) = errSmax_vec;
            [chi,~,~,~,~] = findCritChi3D(U,V,W,desiredAcc);
            chis(k) = max(chi); %extract chi_max at this time and accuracy. 
        end
        
        %Dump data
        if (step <simutimeSteps)
            folder = 'intermediate';
        else
            folder = 'final';
        end
               
        fileName = sprintf('%s/%s/Vortex3D_Re%i_L%.1f_mdt%.0f_mt%.0f.mat'...
            ,outputPath,folder,round(Re),log2(N),1000*dt,1000*t);
        save(fileName,'U','V','W','x','y','z','KX','KY','KZ','Re','dt','t');
                
        %Advance index
        k=k+1;
    end
    
    % Calculate Ekin and vorticity squared (for enstrophy)
    [Ek,~,~, RU, RV, RW] = calcAggQuantities3D(U,V,W);
    Ekin(step+1) = sum(sum(sum(Ek)));
    vortSq(step+1) = sum(sum(sum((RU.^2+RV.^2+RW.^2))))/N^3;
    
    %Perform a single time-step forward
    SolveNavierStokes3D;
    
end

%Output data from the TN analysis
if(lengthScale == round(lengthScale))
    %save(fileName,'U_SingMatMax','U_errT','U_errSmax','V_SingMatMax','V_errT','V_errSmax','W_SingMatMax','W_errT','W_errSmax','Ekin','vortSq','-append');
    save(fileName,'Ekin','vortSq','chis','-append');
else
    save(fileName,'Ekin','vortSq','-append');
end

end