
addpath('../')
addpath('../2DTurbulentFlowProducer')
addpath('../MPSBuilder')

%% Simulation parametres
kmax = 200; %1/h
alpha = 1; %Must be an integer.
T_dds =  0.025;
t_final = 2;
%
mkdir output;
mkdir output/final;
mkdir output/intermediate;
outputPath = 'output/';
%
%%
Re0 = 0.01; log2N = 9; %Re=10.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3; 
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 0.02; log2N = 9; %Re=20.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3; 
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 0.05; log2N = 9; %Re=50.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3; 
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 0.1; log2N = 9; %Re=100.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3; 
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 0.2; log2N = 9; %Re=200.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3; 
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 0.5; log2N = 9; %Re=500.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3; 
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);
%
Re0 = 1; log2N = 10; %Re=1000.
dt = 0.1*2^-(round(log2N)-1);
Re = kmax*Re0*10^3;
simulate2DINSFlow(log2N,kmax,alpha,Re,dt,t_final,T_dds,outputPath);