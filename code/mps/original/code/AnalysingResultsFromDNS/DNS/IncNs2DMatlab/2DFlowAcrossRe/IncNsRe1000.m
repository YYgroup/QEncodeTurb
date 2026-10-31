
addpath('../')
addpath('../2DTurbulentFlowProducer')
addpath('../MPSBuilder')

%% Simulation parametres
kmax = 200; %1/h
alpha = 1; %Must be an integer.
T_dds =  0.025;
t_final = 2.05;
%
mkdir output;
mkdir output/final;
mkdir output/intermediate;
outputPath = 'output/';
%
%%
%
Re0 = 1; log2N = 7; %Re=1000.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3;
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 1; log2N = 8; %Re=1000.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3;
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 1; log2N = log2(362); %Re=1000.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3;
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 1; log2N = 10; %Re=1000.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3;
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);