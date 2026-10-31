function simulate2DINSFlow(lengthScale,kmax,alpha,Re,dt,t_final,TbetweenDDs,outputPath)
%Simulates 2D flow using the scripts in the 2DTurbulentFlowProducer folder
%while simultaneously converting flow fields into a matrix product state type
% tensor network (TN) and extracting and analysing the entanglement spectra.
%t_final, dt and numMakeTN must be integer powers of 10. numMakeTN must be
%less than simutimeStep.

%% Create field
%Define time-evolution variables
N               = round(2^lengthScale); % N*N grid 
simutimeSteps   = ceil(t_final/dt);
StepsBetweenDDs = floor(TbetweenDDs/dt);

% Non-dimensionalise the problem
Umax=1;
L=1;

%Initialise field
CreateFields;

%%
%Initialise time-series arrays
Ekin=zeros(ceil(simutimeSteps/StepsBetweenDDs),1); 
Diss=zeros(ceil(simutimeSteps/StepsBetweenDDs),1); 
maxPhi = zeros(ceil(simutimeSteps/StepsBetweenDDs),1); 
minPhi = zeros(ceil(simutimeSteps/StepsBetweenDDs),1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes simulation begins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if(lengthScale == round(lengthScale))
    U_SingMatMax = zeros(2,floor(simutimeSteps/StepsBetweenDDs));
    V_SingMatMax = zeros(2,floor(simutimeSteps/StepsBetweenDDs));
    %
    D_max = 2^(2*floor(lengthScale/2)); %If even, the max singular value matrix's dimension
    %is 2^lengthScale. If odd, 2^(lengthScale-1).
    U_errT = zeros(D_max,floor(simutimeSteps/StepsBetweenDDs));
    U_errSmax = U_errT;
    %
    V_errT = U_errT;
    V_errSmax = U_errT;
    desiredAcc = [10^-4];
    chis = zeros(1,floor(simutimeSteps/StepsBetweenDDs));
end
k = 1;
%
for step = 0:simutimeSteps
    t = step*dt;
    
    %Data dumps take a long time (approx 4 times as long as performing a
    %time-step of the field forward), so do them sparingly.
    if(mod(step,StepsBetweenDDs)==0 || step == simutimeSteps )
        
        %Represent U,V in TN form and analyse if lengthScale permits
        if(lengthScale == round(lengthScale))
            [SingMats,eEs] = createMPS(U,N);
            [Smax,bond,errSmax_vec,errT_vec] = TNanalyseFlow(SingMats,eEs);
            %
            U_SingMatMax(1,k) = Smax;
            U_SingMatMax(2,k) = bond;
            %
            U_errT(:,k) = errT_vec;
            U_errSmax(:,k) = errSmax_vec;
            
            %Represent V in TN form and analyse
            [SingMats,eEs] = createMPS(V,N);
            [Smax,bond,errSmax_vec,errT_vec] = TNanalyseFlow(SingMats,eEs);
            %
            V_SingMatMax(1,k) = Smax;
            V_SingMatMax(2,k) = bond;
            %
            V_errT(:,k) = errT_vec;
            V_errSmax(:,k) = errSmax_vec;
              [chi,~,~,~,~] = findCritChi2D(U,V,desiredAcc);
              chis(k) = max(chi); %extract chi_max at this time and accuracy. 
        end
        
        %Dump data
        if (step <simutimeSteps)
            folder = 'intermediate';
        else
            folder = 'final';
        end
        
        fileName = sprintf('%s/%s/Vortex2D_Re%i_L%.1f_mdt%.0f_mt%.0f.mat' ...
            ,outputPath,folder,Re,lengthScale,1000*dt,round(1000*t));
        save(fileName,'U','V','Phi','x','y','Re','dt','t');
        
        %Gather data on Energy of flow
        Ekin(k) = sum(sum((U.^2 + V.^2)));%*0.5*dx*dx/L^2;
        maxPhi(k) = max(max(Phi));
        minPhi(k) = min(min(Phi));
        dudx = (real(ifft2( 1i*KX.*fft2(U)))); dudy = (real(ifft2( 1i*KY.*fft2(U)))); 
        dvdx = (real(ifft2( 1i*KX.*fft2(V)))); dvdy = (real(ifft2( 1i*KY.*fft2(V)))); 
        Diss(k) = sum(sum((dudx.^2 + dudy.^2 + dvdx.^2 + dvdy.^2)));%*dx*dx/L^2;
        
        %Advance index
        k=k+1;
    end
    
    %Perform a single time-step forward
    SolveNavierStokes2D;
end

%Output data from the TN analysis
if(lengthScale == round(lengthScale))
    %save(fileName,'U_SingMatMax','U_errT','U_errSmax','V_SingMatMax','V_errT','V_errSmax','Ekin','maxPhi','minPhi','Diss','-append');
    save(fileName,'Ekin','maxPhi','minPhi','Diss','chis','-append');
else
    save(fileName,'Ekin','maxPhi','minPhi','Diss','-append');
end

end