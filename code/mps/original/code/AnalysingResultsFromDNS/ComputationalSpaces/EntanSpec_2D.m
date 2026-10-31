%% First clear and add paths
clear;
addpath('../../plotTools/subplot_tight/subplot_tight');

%% Plot parameters
set(gcf,'color','w');

% Placement
x0= 0.115;
y0= 0.54;
lx = 0.2;
ly = 0.38;

x1=x0-5.5;
y1 = y0+9*50;

%x1=-5;
%y1=+20;

diff_x = lx+0.023;
diff_y = ly+0.015;

figure(1); clf
linWid = 1.5;
mrkSz=4;
fontsz=10;

set(gcf, 'Units','points','Position',  [0, 0, 468, 324+162/3]) % A4 is 595 × 842 points. In arXiv template, default textwidth is 468 pt and text height is 648 pt.

legends_2D = {''};



%% Entanglement Spectrum 2D
load('IncompNS2D/MWS_e2to12_out');

%Pick isocontour line
count = 3; %3 corresponds to 10^-4 L2 squared, <=> 99% accuracy 
%
%% Normalise and clean data
SingValsU250=SingValsU250/sqrt(sum((SingValsU250(1,:)).^2));
SingValsU750=SingValsU750/sqrt(sum((SingValsU750(1,:)).^2));
SingValsU1250=SingValsU1250/sqrt(sum((SingValsU1250(1,:)).^2));
SingValsU1750=SingValsU1750/sqrt(sum((SingValsU1750(1,:)).^2));
SingValsV250=SingValsV250/sqrt(sum((SingValsV250(1,:)).^2));
SingValsV750=SingValsV750/sqrt(sum((SingValsV750(1,:)).^2));
SingValsV1250=SingValsV1250/sqrt(sum((SingValsV1250(1,:)).^2));
SingValsV1750=SingValsV1750/sqrt(sum((SingValsV1750(1,:)).^2));
%
SingValsU250(SingValsU250==0)=NaN;
SingValsU750(SingValsU750==0)=NaN;
SingValsU1250(SingValsU1250==0)=NaN;
SingValsU1750(SingValsU1750==0)=NaN;
SingValsV250(SingValsV250==0)=NaN;
SingValsV750(SingValsV750==0)=NaN;
SingValsV1250(SingValsV1250==0)=NaN;
SingValsV1750(SingValsV1750==0)=NaN;
%
%% Plot Data 2D
xmajticks_2D = 0:2:10;
xminticks_2D = 1:2:9;
xmajticklabels_2D = {'$1$','$4$','$4^2$','$4^3$','$4^4$','$4^5$'};
yticklabels_2D = {'$1$','$2$','$3$','$4$','$5$','$6$','$7$','$8$','$9$'};
contlvls= linspace(log10(10^-17),log10(1),18);

% Plot u_1 @ t=0.25
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0,y0,lx,ly])
subplotEntanSpec(SingValsU250, chi250(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,1,1,1,1);
text(x1,y1,'$u_1$', 'interpreter', 'latex','FontSize',1.2*fontsz);

% Plot u_1 @ t=0.75
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+diff_x,y0,lx,ly])
subplotEntanSpec(SingValsU750, chi750(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,1,0,0,0);

% Plot u_1 @ t=1.25
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+2*diff_x,y0,lx,ly])
subplotEntanSpec(SingValsU1250, chi1250(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,1,0,0,0);

% Plot u_1 @ t=1.75
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+3*diff_x,y0,lx,ly])
subplotEntanSpec(SingValsU1750, chi1750(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,1,0,0,0);

% Plot u_2 @ t=0.25
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV250, chi250(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,0,1,0,0);
text(x1,y1,'$u_2$', 'interpreter', 'latex','FontSize',1.2*fontsz);

% Insert lambda label for colorbar such that it's mostly aligned with the t=0.25 velocity components
lambdaLabel = {'$\frac{\lambda_{\alpha}}{ \sqrt{ \sum_{\alpha} \lambda_{\alpha}^2} }$'};
t = text(x1*1.05,y1*2.33,lambdaLabel);
t.FontSize=1.2*fontsz;
t.Interpreter='latex';

% Plot u_2 @ t=0.75
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+diff_x,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV750, chi750(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,0,0,0,0);

% Plot u_2 @ t=1.25
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+2*diff_x,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV1250, chi1250(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,0,0,0,0);

% Plot u_2 @ t=1.75
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+3*diff_x,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV1750, chi1750(:,count), contlvls, fontsz,linWid,xmajticks_2D,xmajticklabels_2D,xminticks_2D,yticklabels_2D,0,0,0,0);

%% Colorscale and colorbar
colormap(jet(length(contlvls)-1));
caxis([log10(10^-17) log10(1)]);
ticklabels = {'','$10^{-16}$','','$10^{-14}$','','$10^{-12}$','','$10^{-10}$','','$10^{-8}$','','$10^{-6}$','','$10^{-4}$','','$10^{-2}$','','$1$'};
h=colorbar('Location','southoutside','Ticks',contlvls,'TickLabels',ticklabels,'FontSize',fontsz);
set( h, 'YDir', 'reverse' );
set(h, 'Position', [x0 .05 .87 .07])
h.LineWidth=linWid;
h.TickLabelInterpreter = 'latex';

%% Save
print(gcf,'../../../results/EntanSpec_2D.png','-dpng','-r450')
%print -depsc2 EntanSpec_2D