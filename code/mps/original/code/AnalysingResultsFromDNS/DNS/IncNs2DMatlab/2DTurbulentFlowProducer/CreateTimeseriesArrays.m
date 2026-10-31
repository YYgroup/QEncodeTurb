%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example to create arrays for saving 
% kinetic energy and 
% dissipation rate time series data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ekin=zeros(simutimeSteps,1); 
Diss=zeros(simutimeSteps,1); 
Time=linspace(0,simutimeSeconds,simutimeSteps)';
