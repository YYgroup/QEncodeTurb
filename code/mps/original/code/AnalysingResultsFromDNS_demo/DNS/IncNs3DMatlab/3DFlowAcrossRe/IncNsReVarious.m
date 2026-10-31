
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
Re = 2*pi*10; %Re=10.
log2N = 4;
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
Re = 2*pi*20; %Re=20.
log2N = 5;
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
Re = 2*pi*50; %Re=50.
log2N = 6;
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
Re = 2*pi*100; %Re=100.
log2N = 6;
dt = 0.1*2^-(ceil(log2N)-2);
simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
%
% Re = 2*pi*200; %Re=200.
% log2N = 7;
% dt = 0.1*2^-(ceil(log2N)-2);
% simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
% %
% Re = 2*pi*500; %Re=500.
% log2N = 8;
% dt = 0.1*2^-(ceil(log2N)-2);
% simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);
% %
% Re = 2*pi*800; %Re=800.
% log2N = 8;
% dt = 0.1*2^-(ceil(log2N)-2);
% simulate3DINSFlow(log2N,kmax,Re,dt,t_final,T_dds,outputPath);