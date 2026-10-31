%% First clear and add paths
clear;

%% chi_1% for 2D & 3D.
%
load('../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/final/Vortex2D_Re2000_L5.0_mdt6_mt2000.mat'); %Re=10.
chi2D_Re10 = max(chis);
load('../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/final/Vortex2D_Re4000_L6.0_mdt3_mt2000.mat'); %Re=20.
chi2D_Re20 = max(chis);
load('../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/final/Vortex2D_Re10000_L7.0_mdt2_mt2000.mat'); %Re=50.
chi2D_Re50 = max(chis);
load('../DNS/IncNs2DMatlab/2DFlowAcrossRe/output/final/Vortex2D_Re20000_L8.0_mdt1_mt2000.mat'); %Re=100.
chi2D_Re100 = max(chis);
%
load('../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/final/Vortex3D_Re63_L4.0_mdt25_mt2000.mat'); %Re=10.
chi3D_Re10 = max(chis);
load('../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/final/Vortex3D_Re126_L5.0_mdt12_mt2000.mat'); %Re=20.
chi3D_Re20 = max(chis);
load('../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/final/Vortex3D_Re314_L6.0_mdt6_mt2000.mat'); %Re=50.
chi3D_Re50 = max(chis);
load('../DNS/IncNs3DMatlab/3DFlowAcrossRe/output/final/Vortex3D_Re628_L6.0_mdt6_mt2000.mat'); %Re=100.
chi3D_Re100 = max(chis);
%
chis2D = [chi2D_Re10,chi2D_Re20,chi2D_Re50,chi2D_Re100];
chis3D = [chi3D_Re10,chi3D_Re20,chi3D_Re50,chi3D_Re100];
%
%% Save
Re2D = ([10,20,50,100]);
Re3D = ([10,20,50,100]);

chis2D_log2 = log2(chis2D);%/log2(8);
chis3D_log2 = log2(chis3D);%/log2(8);

%Estimate exponent of TGV power-law
X = ones(4,2); X(:,2) = log10(Re3D');
Y = chis3D_log2';
B = X\Y;
b0 = B(1);
b1 = B(2);

nikEstimate3D_log2 = b0 + b1*log10(Re3D);


save('IncompNS_ChisRe.mat','nikEstimate3D_log2','chis2D_log2','chis3D_log2','Re2D','Re3D');
