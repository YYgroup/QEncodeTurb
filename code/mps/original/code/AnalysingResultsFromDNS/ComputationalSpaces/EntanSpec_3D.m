%% First clear and add paths
clear;
addpath('../../plotTools/subplot_tight/subplot_tight');

%% Plot parameters
set(gcf,'color','w');

% Placement
x0= 0.115;
y0= 0.674;
lx = 0.202;
ly = 0.415/1.5;

diff_x = lx+0.021;
diff_y = ly+0.01;

x1 = x0-6.7;
y1 = y0+8*50;

%x1 = -6;
%y1=+17;

figure(1); clf
linWid = 1.5;
mrkSz=4;
fontsz=10;

set(gcf, 'Units','points','Position',  [0, 0, 468, 324*1.5+162/3]) % A4 is 595 × 842 points. In arXiv template, default textwidth is 468 pt and text height is 648 pt.

%% Entanglement Spectrum 3D
load('IncompNS3D/MWS_e2to12_out');

%Pick isocontour line
count = 3; %3 corresponds to 10^-4 L2 squared, <=> 99% accuracy 
%
%% Normalise and clean data
SingValsU200=SingValsU200/sqrt(sum((SingValsU200(1,:)).^2));
SingValsU800=SingValsU800/sqrt(sum((SingValsU800(1,:)).^2));
SingValsU1400=SingValsU1400/sqrt(sum((SingValsU1400(1,:)).^2));
SingValsU2000=SingValsU2000/sqrt(sum((SingValsU2000(1,:)).^2));
SingValsV200=SingValsV200/sqrt(sum((SingValsV200(1,:)).^2));
SingValsV800=SingValsV800/sqrt(sum((SingValsV800(1,:)).^2));
SingValsV1400=SingValsV1400/sqrt(sum((SingValsV1400(1,:)).^2));
SingValsV2000=SingValsV2000/sqrt(sum((SingValsV2000(1,:)).^2));
SingValsW200=SingValsW200/sqrt(sum((SingValsW200(1,:)).^2));
SingValsW800=SingValsW800/sqrt(sum((SingValsW800(1,:)).^2));
SingValsW1400=SingValsW1400/sqrt(sum((SingValsW1400(1,:)).^2));
SingValsW2000=SingValsW2000/sqrt(sum((SingValsW2000(1,:)).^2));
%
SingValsU200(SingValsU200==0)=NaN;
SingValsU800(SingValsU800==0)=NaN;
SingValsU1400(SingValsU1400==0)=NaN;
SingValsU2000(SingValsU2000==0)=NaN;
SingValsV200(SingValsV200==0)=NaN;
SingValsV800(SingValsV800==0)=NaN;
SingValsV1400(SingValsV1400==0)=NaN;
SingValsV2000(SingValsV2000==0)=NaN;
SingValsW200(SingValsW200==0)=NaN;
SingValsW800(SingValsW800==0)=NaN;
SingValsW1400(SingValsW1400==0)=NaN;
SingValsW2000(SingValsW2000==0)=NaN;
%
%% Plot Data 3D
xmajticks_3D = 0:3:12;
xminticks_3D = [1,2,4,5,7,8,10,11];
xmajticklabels_3D = {'$1$','$8$','$8^2$','$8^3$','$8^4$'};
yticklabels_3D = {'$1$','$2$','$3$','$4$','$5$','$6$','$7$'};
contlvls= linspace(log10(10^-17),log10(1),18);

% Plot u_1 @ t=0.2
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0,y0,lx,ly])
subplotEntanSpec(SingValsU200, chi200(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,1,1,1,1);
text(x1,y1,'$u_1$', 'interpreter', 'latex','FontSize',1.2*fontsz);

% Plot u_1 @ t=0.8
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+diff_x,y0,lx,ly])
subplotEntanSpec(SingValsU800, chi800(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,1,0,0,0);


% Plot u_1 @ t=1.4
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+2*diff_x,y0,lx,ly])
subplotEntanSpec(SingValsU1400, chi1400(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,1,0,0,0);

% Plot u_1 @ t=2
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+3*diff_x,y0,lx,ly])
subplotEntanSpec(SingValsU2000, chi2000(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,1,0,0,0);

% Plot u_2 @ t=0.2
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV200, chi200(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,1,0,0);
text(x1,y1,'$u_2$', 'interpreter', 'latex','FontSize',1.2*fontsz);

% Plot u_2 @ t=0.8
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+diff_x,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV800,chi800(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,0,0,0);

% Plot u_2 @ t=1.4
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+2*diff_x,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV1400,chi1400(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,0,0,0);

% Plot u_2 @ t=2
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+3*diff_x,y0-diff_y,lx,ly])
subplotEntanSpec(SingValsV2000, chi2000(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,0,0,0);

% Plot u_3 @ t=0.2
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0,y0-diff_y*2,lx,ly])
subplotEntanSpec(SingValsW200, chi200(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,1,0,0);
text(x1,y1,'$u_3$', 'interpreter', 'latex','FontSize',1.2*fontsz);

% Insert lambda label for colorbar such that it's mostly aligned with the t=0.2 velocity components
lambdaLabel = {'$\frac{\lambda_{\alpha}}{ \sqrt{ \sum_{\alpha} \lambda_{\alpha}^2}}$'};
t = text(x1*1.02,y1*2.1,lambdaLabel);
t.FontSize=1.2*fontsz;
t.Interpreter='latex';

% Plot u_3 @ t=0.8
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+diff_x,y0-diff_y*2,lx,ly])
subplotEntanSpec(SingValsW800, chi800(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,0,0,0);

% Plot u_3 @ t=1.4
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+2*diff_x,y0-diff_y*2,lx,ly])
subplotEntanSpec(SingValsW1400, chi1400(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,0,0,0);

% Plot u_3 @ t=2
subplot_tight(5,3,1,[0.00 0.00],'Position', [x0+3*diff_x,y0-diff_y*2,lx,ly])
subplotEntanSpec(SingValsW2000, chi2000(:,count), contlvls, fontsz,linWid,xmajticks_3D,xmajticklabels_3D,xminticks_3D,yticklabels_3D,0,0,0,0);


%% Colorscale and colorbar
colormap(jet(length(contlvls)-1));
caxis([log10(10^-17) log10(1)]);
ticklabels = {'','$10^{-16}$','','$10^{-14}$','','$10^{-12}$','','$10^{-10}$','','$10^{-8}$','','$10^{-6}$','','$10^{-4}$','','$10^{-2}$','','$1$'};
h=colorbar('Location','southoutside','Ticks',contlvls,'TickLabels',ticklabels,'FontSize',fontsz);
set(h, 'Position', [.86 .01 .07 (3*ly+0.02)])
set( h, 'YDir', 'reverse' );
set(h, 'Position', [x0 .03 .87 .07/1.5])
h.LineWidth=linWid;
h.TickLabelInterpreter = 'latex';

%% Save
print(gcf,'../../../results/EntanSpec_3D.png','-dpng','-r450')
%print -depsc2 EntanSpec_3D