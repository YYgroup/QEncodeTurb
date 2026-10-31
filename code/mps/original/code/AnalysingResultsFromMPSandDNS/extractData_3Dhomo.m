%% Compares simulation results of TNT and Matlab 

%% First add paths to access required functions and files
clear all; warning('off','all');

addpath('../PDE_matfiles/General');
addpath('../PDE_matfiles/MPSUtility');
addpath('../PDE_matfiles/MPSBuilderNik');
addpath('../PDE_matfiles/NS3DIncomp');

TNTpath_incompNS = '../../data/IncompNs3D/';
MATpath_incompNS = '../AnalysingResultsFromDNS/DNS/IncNs3DMatlab/3DFlowAcrossRe/output/';


%% dedt & enst data
% For Article
filenameStats_L8C192Re800 = strcat(TNTpath_incompNS,sprintf('NsIncomp3DdedtenstL8C192Re5027.mat'));
filenameStats_L8C128Re800 = strcat(TNTpath_incompNS,sprintf('NsIncomp3DEkinEnstL8C128Re5027.mat'));
filenameStats_L8C96Re800  = strcat(TNTpath_incompNS,sprintf('NsIncomp3DEkinEnstL8C96Re5027.mat'));

%% data for visualisation

%visualisation:
%DNS
filenameFieldt200_L8MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L8.0_mdt2_mt200.mat'));
filenameFieldt800_L8MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L8.0_mdt2_mt800.mat'));
filenameFieldt1400_L8MRe800 =strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L8.0_mdt2_mt1400.mat'));
filenameFieldt2000_L8MRe800 =strcat(MATpath_incompNS,sprintf('/final/Vortex3D_Re5027_L8.0_mdt2_mt2000.mat'));

%88^3
filenameFieldt200_N88MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L6.5_mdt3_mt200.mat'));
filenameFieldt800_N88MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L6.5_mdt3_mt800.mat'));
filenameFieldt1400_N88MRe800 =  strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L6.5_mdt3_mt1400.mat'));
filenameFieldt2000_N88MRe800 = strcat(MATpath_incompNS,sprintf('/final/Vortex3D_Re5027_L6.5_mdt3_mt2000.mat'));

%70^3
filenameFieldt200_N70MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L6.1_mdt3_mt200.mat'));
filenameFieldt800_N70MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L6.1_mdt3_mt800.mat'));
filenameFieldt1400_N70MRe800 =  strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L6.1_mdt3_mt1400.mat'));
filenameFieldt2000_N70MRe800 = strcat(MATpath_incompNS,sprintf('/final/Vortex3D_Re5027_L6.1_mdt3_mt2000.mat'));

%60^3
filenameFieldt200_N60MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L5.9_mdt6_mt200.mat'));
filenameFieldt800_N60MRe800 = strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L5.9_mdt6_mt800.mat'));
filenameFieldt1400_N60MRe800 =  strcat(MATpath_incompNS,sprintf('/intermediate/Vortex3D_Re5027_L5.9_mdt6_mt1400.mat'));
filenameFieldt2000_N60MRe800 = strcat(MATpath_incompNS,sprintf('/final/Vortex3D_Re5027_L5.9_mdt6_mt2000.mat'));

%chi=192
filenameFieldt200_L8C192Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C192Re5027mt200.mat'));
filenameFieldt800_L8C192Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C192Re5027mt800.mat'));
filenameFieldt1400_L8C192Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C192Re5027mt1400.mat'));
filenameFieldt2000_L8C192Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C192Re5027mt2000.mat'));

%chi=128
filenameFieldt200_L8C128Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C128Re5027mt200.mat'));
filenameFieldt800_L8C128Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C128Re5027mt800.mat'));
filenameFieldt1400_L8C128Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C128Re5027mt1400.mat'));
filenameFieldt2000_L8C128Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C128Re5027mt2000.mat'));

%chi=96
filenameFieldt200_L8C96Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C96Re5027mt200.mat'));
filenameFieldt800_L8C96Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C96Re5027mt800.mat'));
filenameFieldt1400_L8C96Re800 = strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C96Re5027mt1400.mat'));
filenameFieldt2000_L8C96Re800 =strcat(TNTpath_incompNS,sprintf('NsIncomp3DFieldOutL8C96Re5027mt2000.mat'));


%% Load dedt&enst data for article
%Loading of C192 energy decay and enstropy data is special because this simulation crashed
%and had to be restarted, and thus I had to reconnect the dedt&denst arrays
%manually. I then took the opportunity to also manually process it (i.e.
%calculate derivatives and scale enstropy with Reynolds number).
if(fopen(filenameStats_L8C192Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameStats_L8C192Re800);
else
    load(filenameStats_L8C192Re800);
end
%C128
if(fopen(filenameStats_L8C128Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameStats_L8C128Re800);
else
    load(filenameStats_L8C128Re800);
    dedt_L8C128Re800 = diff(Ekin)/dt;
    enst_L8C128Re800 = (Enst(1:end-1)+Enst(2:end))/Re;
    dt_L8C128Re800=dt;
end
%C96
if(fopen(filenameStats_L8C96Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameStats_L8C96Re800);
else
    load(filenameStats_L8C96Re800);
    dedt_L8C96Re800 = diff(Ekin)/dt;
    enst_L8C96Re800 = (Enst(1:end-1)+Enst(2:end))/Re;
    dt_L8C96Re800=dt;
end
%N256
if(fopen(filenameFieldt2000_L8MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_L8MRe800);
else
    load(filenameFieldt2000_L8MRe800);
    dedt_L8MRe800 = diff(Ekin)/dt;
    enst_L8MRe800 = (vortSq(1:end-1)+vortSq(2:end))/Re;
    dt_L8MRe800=dt;
end
%N88
if(fopen(filenameFieldt2000_N88MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_N88MRe800);
else
    load(filenameFieldt2000_N88MRe800);
    dedt_N88MRe800 = diff(Ekin)/dt;
    enst_N88MRe800 = (vortSq(1:end-1)+vortSq(2:end))/Re;
    dt_N88MRe800=dt;
end
%N70
if(fopen(filenameFieldt2000_N70MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_N70MRe800);
else
    load(filenameFieldt2000_N70MRe800);
    dedt_N70MRe800 = diff(Ekin)/dt;
    enst_N70MRe800 = (vortSq(1:end-1)+vortSq(2:end))/Re;
    dt_N70MRe800=dt;
end
%N60
if(fopen(filenameFieldt2000_N60MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_N60MRe800);
else
    load(filenameFieldt2000_N60MRe800);
    dedt_N60MRe800 = diff(Ekin)/dt;
    enst_N60MRe800 = (vortSq(1:end-1)+vortSq(2:end))/Re;
    dt_N60MRe800=dt;
end

%% Load DNS and URDNS field data for visualisation

%N256, i.e. DNS
if(fopen(filenameFieldt200_L8MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_L8MRe800);
else
    load(filenameFieldt200_L8MRe800);
    lambda2_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_L8MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_L8MRe800);
else
    load(filenameFieldt800_L8MRe800);
    lambda2_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_L8MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_L8MRe800);
else
    load(filenameFieldt1400_L8MRe800);
    lambda2_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_L8MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_L8MRe800);
else
    load(filenameFieldt2000_L8MRe800);
    lambda2_t2000 = getLambda2(U, V, W);
end
%
%N88
if(fopen(filenameFieldt200_N88MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_N88MRe800);
else
    load(filenameFieldt200_N88MRe800);
    lambda2N88_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_N88MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_N88MRe800);
else
    load(filenameFieldt800_N88MRe800);
    lambda2N88_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_N88MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_N88MRe800);
else
    load(filenameFieldt1400_N88MRe800);
    lambda2N88_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_N88MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_N88MRe800);
else
    load(filenameFieldt2000_N88MRe800);
    lambda2N88_t2000 = getLambda2(U, V, W);
end
%
%N70
if(fopen(filenameFieldt200_N70MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_N70MRe800);
else
    load(filenameFieldt200_N70MRe800);
    lambda2N70_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_N70MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_N70MRe800);
else
    load(filenameFieldt800_N70MRe800);
    lambda2N70_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_N70MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_N70MRe800);
else
    load(filenameFieldt1400_N70MRe800);
    lambda2N70_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_N70MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_N70MRe800);
else
    load(filenameFieldt2000_N70MRe800);
    lambda2N70_t2000 = getLambda2(U, V, W);
end
%
%N60
if(fopen(filenameFieldt200_N60MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_N60MRe800);
else
    load(filenameFieldt200_N60MRe800);
    lambda2N60_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_N60MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_N60MRe800);
else
    load(filenameFieldt800_N60MRe800);
    lambda2N60_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_N60MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_N60MRe800);
else
    load(filenameFieldt1400_N60MRe800);
    lambda2N60_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_N60MRe800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_N60MRe800);
else
    load(filenameFieldt2000_N60MRe800);
    lambda2N60_t2000 = getLambda2(U, V, W);
end
%

%% Load MPS field data for visualisation

%chi=192
if(fopen(filenameFieldt200_L8C192Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_L8C192Re800);
else
    load(filenameFieldt200_L8C192Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C192_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_L8C192Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_L8C192Re800);
else
    load(filenameFieldt800_L8C192Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C192_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_L8C192Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_L8C192Re800);
else
    load(filenameFieldt1400_L8C192Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C192_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_L8C192Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_L8C192Re800);
else
    load(filenameFieldt2000_L8C192Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C192_t2000 = getLambda2(U, V, W);
end
%
%chi=128
if(fopen(filenameFieldt200_L8C128Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_L8C128Re800);
else
    load(filenameFieldt200_L8C128Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C128_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_L8C128Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_L8C128Re800);
else
    load(filenameFieldt800_L8C128Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C128_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_L8C128Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_L8C128Re800);
else
    load(filenameFieldt1400_L8C128Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C128_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_L8C128Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_L8C128Re800);
else
    load(filenameFieldt2000_L8C128Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C128_t2000 = getLambda2(U, V, W);
end
%
%chi=96
if(fopen(filenameFieldt200_L8C96Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt200_L8C96Re800);
else
    load(filenameFieldt200_L8C96Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C96_t200 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt800_L8C96Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt800_L8C96Re800);
else
    load(filenameFieldt800_L8C96Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C96_t800 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt1400_L8C96Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt1400_L8C96Re800);
else
    load(filenameFieldt1400_L8C96Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C96_t1400 = getLambda2(U, V, W);
end
%
if(fopen(filenameFieldt2000_L8C96Re800) == -1)
    fprintf('Error: file "%s" doesn''t exist.\n',filenameFieldt2000_L8C96Re800);
else
    load(filenameFieldt2000_L8C96Re800);
    [U,V,W]= reconstructFields3D(u1,u2,u3);
    lambda2C96_t2000 = getLambda2(U, V, W);
end

%% clear unneeded data and save
clear U V W KX KY KZ;
%
save('extractedData_3Dhomo');

