
addpath('../')
addpath('../3DTurbulentFlowProducer')
addpath('../MPSBuilder')

%% Simulation Parameters 
kmax = 2*pi;
T_dds = 0.025;
t_final =  2;%0.1;
%
mkdir output;
mkdir output/final;
mkdir output/intermediate;
outputPath = 'output/';

%%

Re = 2*pi*800; %Re=800.
log2N = log2(60);
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
Re = 2*pi*800; %Re=800.
log2N = log2(70);
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
Re = 2*pi*800; %Re=800.
log2N = log2(88);
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
Re = 2*pi*800; %Re=800.
log2N = 8;
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);