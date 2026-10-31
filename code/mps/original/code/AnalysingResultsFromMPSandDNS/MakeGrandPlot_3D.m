%clear; clc;
restoredefaultpath;
clear restoredefaultpath;

addpath('../PDE_matfiles/MPSUtility');
addpath('../PDE_matfiles/MPSBuilderNik');
addpath('../PDE_matfiles/NS3DIncomp');

addpath('../plotTools/subplot_tight/subplot_tight');
addpath('../plotTools/legendflex');

set(0,'DefaultLegendAutoUpdate','off')

load('extractedData_3Dhomo.mat')

set(gcf,'color','w');

% Placement
x0= 0.155;
y0= 0.18;
X0 = 0.6;
Y0 = y0;

lx = 0.65;
ly = 0.8;
Lx = lx;
Ly = ly;

xLgd=0.13;
yLgd=0.59;

lgd_width=[xLgd yLgd 0.2 0.2];
lgd_pos=[-75 -14];

fontSz = 12;


figure(1); clf
legends = {'$\mathrm{DNS}$','$\chi{=}192$','$\chi{=}128$','$\chi{=}96$','$88^3\textnormal{ grid}$','$70^3\textnormal{ grid}$','$60^3\textnormal{ grid}$'};
linWid = 1.5;
linWidAx= 1.5;
mrkSz=7;

set(gcf, 'Units','points','Position',  [2500, 500, 330, 200]) % A4 is 595 × 842 points. In arXiv template, default textwidth is 468 pt and text height is 648 pt.



%% Enstrophy
%
mainFig=figure(1);
subplot_tight(1,2,1,[0.00 0.00],'Position',[x0,y0,lx,ly])

Ns = length(dedt_L8C96Re800)+1;
Tf = dt_L8C96Re800*(Ns-1);
t=linspace(dt_L8C96Re800/2,Tf-dt_L8C96Re800/2,Ns-1);
Ns_88_80_70 = length(dedt_N88MRe800)+1;
Tf_88_80_70 = dt_N88MRe800*(Ns_88_80_70-1);
t_88_80_70=linspace(dt_N88MRe800/2,Tf_88_80_70-dt_N88MRe800/2,Ns_88_80_70-1);
Ns_60 = length(dedt_N60MRe800)+1;
Tf_60 = dt_N60MRe800*(Ns_60-1);
t_60=linspace(dt_N60MRe800/2,Tf_60-dt_N88MRe800/2,Ns_60-1);

currPlot(1)=plot(t,-0.5*dedt_L8MRe800/Ekin(1),'LineWidth',1.5*linWid,'color', 'black'); hold on;
currPlot(2)=plot(t,-dedt_L8C192Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0,1]);
currPlot(3)=plot(t,-dedt_L8C128Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0.5,0]);
currPlot(4)=plot(t,-dedt_L8C96Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [1,0,0]);
currPlot(5)=plot(t_88_80_70,-.5*dedt_N88MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0,0.75,.75]);
currPlot(6)=plot(t_88_80_70,-.5*dedt_N70MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0.75,0,0.75]);
currPlot(7)=plot(t_60,-.5*dedt_N60MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[.75,.75,0]);
%
currPlot2(1)=plot(t(1:2^6:end),+0.5*enst_L8MRe800(1:2^6:end)/Ekin(1),'*','color', 'black','MarkerSize',mrkSz);
currPlot2(2)=plot(t(1:2^6:end),+enst_L8C192Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0,1],'MarkerSize',mrkSz);
currPlot2(3)=plot(t(1:2^6:end),+enst_L8C128Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0.5,0],'MarkerSize',mrkSz);
currPlot2(4)=plot(t(1:2^6:end),+enst_L8C96Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [1,0,0],'MarkerSize',mrkSz);
currPlot2(5)=plot(t_88_80_70(1:2^6:end),+.5*enst_N88MRe800(1:2^6:end)/Ekin(1),'.','color',[0,0.75,.75],'MarkerSize',1.5*mrkSz); hold on;
currPlot2(6)=plot(t_88_80_70(1:2^6:end),+.5*enst_N70MRe800(1:2^6:end)/Ekin(1),'.','color',[0.75,0,0.75],'MarkerSize',1.5*mrkSz); hold on;
currPlot2(7)=plot(t_60(1:2^5:end),+.5*enst_N60MRe800(1:2^5:end)/Ekin(1),'.','color',[.75,.75,0],'MarkerSize',1.5*mrkSz); hold on;

xticks([0 .125 .25 .375 .5 .625 .75 .875 1 1.125 1.25 1.375 1.5 1.625 1.75 1.875 2] );yticks([0,0.05,0.1,0.15,0.2,.25,0.3,.35,0.4,.45,0.5,.55,0.6,.65,.7] ); axis([0 2 0 .675])
xticklabels({'$0$','','$0.25$','','$0.5$','','$0.75$','','$1$','','$1.25$','','$1.5$','','$1.75$','','$2$'});
yticklabels({'$0$','','$0.1$','','$0.2$','','$0.3$','','$0.4$','','$0.5$','','$0.6$','','$0.7$'});

set(gca,'TickLabelInterpreter', 'latex');
xlabel('$t/T_0$','interpreter', 'latex');
ax=gca;
ax.LineWidth=linWidAx; set(gca,'layer','top');
ax.FontSize = fontSz;
ylabel('$\frac{\varepsilon(t)}{E_0/T_0}$, $\frac{\zeta(t)}{E_0/T_0}$','interpreter','latex','fontsize',1.5*fontSz);

[lgd,h_legend] = legendflex(currPlot,legends, 'xscale',0.4,'ref', gcf, 'anchor', [2,2], 'buffer',lgd_pos, 'nrow',7,'ncol',1, 'FontSize',fontSz,'Position', lgd_width,'interpreter','latex','Box','off');
grid off;

chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)])
%
%% Inline Figures
rec1=rectangle('position',[1.2 0.51 .35 .151],'LineWidth',linWidAx);
rec2=rectangle('position',[1.88 0.3 .1 .06],'LineWidth',linWidAx);
rec3=rectangle('position',[1.6 0.45 .15 .08],'LineWidth',linWidAx);

rec1.EdgeColor(4) = 0.7; 
rec2.EdgeColor(4) = 0.7;
rec3.EdgeColor(4) = 0.7;

%
linje1=line([1.33 1.2], [.51 .25],'linewidth',linWidAx,'color','black','LineStyle',':');
linje1.Color(4) = 0.7;
linje2=annotation('line',[.8 0.84], [0.58 .5],'LineWidth',linWidAx,'LineStyle',':');
linje2.Color = [0.3,0.3,0.3];
linje3=annotation('line',[.723 0.84], [0.75 .85],'LineWidth',linWidAx,'LineStyle',':');
linje3.Color = [0.3,0.3,0.3];

%
%create axes 
ax1=axes('Parent',mainFig,'Position',[x0+0.27 y0+0.02 .33 0.36]);
ax2=axes('Parent',mainFig,'Position',[x0+0.69 y0-0.0 .15 0.4]);
ax3=axes('Parent',mainFig,'Position',[x0+0.69 y0+0.45 .15 0.35]);
%

%% Inline 1
box(ax1); hold(ax1);

plot(ax1,t,-0.5*dedt_L8MRe800/Ekin(1),'LineWidth',1*linWid,'color', 'black');
currPlot=plot(ax1,t,-dedt_L8C192Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0,1]);
% currPlot.Color(4)=.85;
currPlot=plot(ax1,t,-dedt_L8C128Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0.5,0]);
%currPlot.Color(4)=.7;
currPlot=plot(ax1,t,-dedt_L8C96Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [1,0,0]);
%currPlot.Color(4)=.5;
currPlot=plot(ax1,t_88_80_70,-.5*dedt_N88MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0,0.75,.75]);
%currPlot.Color(4)=.85;
currPlot=plot(ax1,t_88_80_70,-.5*dedt_N70MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0.75,0,0.75]);
%currPlot.Color(4)=.7;
currPlot=plot(ax1,t_60,-.5*dedt_N60MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[.75,.75,0]);
%currPlot.Color(4)=.5;
%
plot(ax1,t(1:2^6:end),+0.5*enst_L8MRe800(1:2^6:end)/Ekin(1),'*','color', 'black','MarkerSize',mrkSz);
plot(ax1,t(1:2^6:end),+enst_L8C192Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0,1],'MarkerSize',mrkSz);
plot(ax1,t(1:2^6:end),+enst_L8C128Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0.5,0],'MarkerSize',mrkSz);
plot(ax1,t(1:2^6:end),+enst_L8C96Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [1,0,0],'MarkerSize',mrkSz);
plot(ax1,t_88_80_70(1:2^5:end),+.5*enst_N88MRe800(1:2^5:end)/Ekin(1),'.','color',[0,0.75,.75],'MarkerSize',1.5*mrkSz);
plot(ax1,t_88_80_70(1:2^5:end),+.5*enst_N70MRe800(1:2^5:end)/Ekin(1),'.','color',[0.75,0,0.75],'MarkerSize',1.5*mrkSz);
plot(ax1,t_60(1:2^4:end),+.5*enst_N60MRe800(1:2^4:end)/Ekin(1),'.','color',[.75,.75,0],'MarkerSize',1.5*mrkSz);
%
axis(ax1,[1.2 1.55 0.51 0.661]); 
%grid(ax1);
xticks(ax1,[1.25 1.375 1.5 1.625]); yticks(ax1, [.55 .6 .65]);
xticklabels(ax1,{'' ''}); yticklabels(ax1,{'' '' ''})
ax1.LineWidth = linWidAx;
%
grid off;

%% Inline 2
box(ax2); hold(ax2);
plot(ax2,t,-0.5*dedt_L8MRe800/Ekin(1),'LineWidth',linWid,'color', 'black');
plot(ax2,t(1:2^6:end),+0.5*enst_L8MRe800(1:2^6:end)/Ekin(1),'*','color', 'black','MarkerSize',mrkSz);
currPlot=plot(ax2,t,-dedt_L8C192Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0,1]);
%currPlot.Color(4)=.85;
plot(ax2,t(1:2^6:end),+enst_L8C192Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0,1],'MarkerSize',mrkSz);
currPlot=plot(ax2,t,-dedt_L8C128Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0.5,0]);
%currPlot.Color(4)=.7;
plot(ax2,t(1:2^6:end),+enst_L8C128Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0.5,0],'MarkerSize',mrkSz);
currPlot=plot(ax2,t,-dedt_L8C96Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [1,0,0]);
%currPlot.Color(4)=.5;
plot(ax2,t(1:2^6:end),+enst_L8C96Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [1,0,0],'MarkerSize',mrkSz);
currPlot=plot(ax2,t_88_80_70,-.5*dedt_N88MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0,0.75,.75]);
%currPlot.Color(4)=.85;
plot(ax2,t_88_80_70(1:2^5:end),+.5*enst_N88MRe800(1:2^5:end)/Ekin(1),'.','color',[0,0.75,.75],'MarkerSize',1.5*mrkSz);
currPlot=plot(ax2,t_88_80_70,-.5*dedt_N70MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0.75,0,0.75]);
%currPlot.Color(4)=.7;
plot(ax2,t_88_80_70(1:2^5:end),+.5*enst_N70MRe800(1:2^5:end)/Ekin(1),'.','color',[0.75,0,0.75],'MarkerSize',1.5*mrkSz);
currPlot=plot(ax2,t_60,-.5*dedt_N60MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[.75,.75,0]);
%currPlot.Color(4)=.5;
plot(ax2,t_60(1:2^4:end),+.5*enst_N60MRe800(1:2^4:end)/Ekin(1),'.','color',[.75,.75,0],'MarkerSize',1.5*mrkSz);
%
axis(ax2,[1.88 1.98 0.3 0.36]); 
%grid(ax2);
xticks(ax2,[1.75 1.875 2]); yticks(ax2, [.25 .3 .35 .4]);
xticklabels(ax2,{'' ''}); yticklabels(ax2,{'' '' ''})
ax2.LineWidth = linWidAx;
%
%% Inline 3
box(ax3); hold(ax3);
plot(ax3,t,-0.5*dedt_L8MRe800/Ekin(1),'LineWidth',linWid,'color', 'black');
plot(ax3,t(1:2^6:end),+0.5*enst_L8MRe800(1:2^6:end)/Ekin(1),'*','color', 'black','MarkerSize',mrkSz);
currPlot=plot(ax3,t,-dedt_L8C192Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0,1]);
%currPlot.Color(4)=.85;
plot(ax3,t(1:2^6:end),+enst_L8C192Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0,1],'MarkerSize',mrkSz);
currPlot=plot(ax3,t,-dedt_L8C128Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [0,0.5,0]);
%currPlot.Color(4)=.7;
plot(ax3,t(1:2^6:end),+enst_L8C128Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [0,0.5,0],'MarkerSize',mrkSz);
currPlot=plot(ax3,t,-dedt_L8C96Re800(1:Ns-1)/Ekin(1),'LineWidth',linWid,'color', [1,0,0]);
%currPlot.Color(4)=.5;
plot(ax3,t(1:2^6:end),+enst_L8C96Re800(1:2^6:Ns-1)/Ekin(1),'x','color', [1,0,0],'MarkerSize',mrkSz);
currPlot=plot(ax3,t_88_80_70,-.5*dedt_N88MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0,0.75,.75]);
%currPlot.Color(4)=.85;
plot(ax3,t_88_80_70(1:2^5:end),+.5*enst_N88MRe800(1:2^5:end)/Ekin(1),'.','color',[0,0.75,.75],'MarkerSize',1.5*mrkSz);
currPlot=plot(ax3,t_88_80_70,-.5*dedt_N70MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[0.75,0,0.75]);
%currPlot.Color(4)=.7;
plot(ax3,t_88_80_70(1:2^5:end),+.5*enst_N70MRe800(1:2^5:end)/Ekin(1),'.','color',[0.75,0,0.75],'MarkerSize',1.5*mrkSz);
currPlot=plot(ax3,t_60,-.5*dedt_N60MRe800/Ekin(1),'LineWidth',linWid,'LineStyle','-.','color',[.75,.75,0]);
%currPlot.Color(4)=.5;
plot(ax3,t_60(1:2^4:end),+.5*enst_N60MRe800(1:2^4:end)/Ekin(1),'.','color',[.75,.75,0],'MarkerSize',1.5*mrkSz);
%
axis(ax3,[1.6 1.75 0.45 0.53]); 
%grid(ax3);
xticks(ax3,[1.5 1.625 1.75 1.875 2]); yticks(ax3, [.25 .3 .35 .4 .45 .5]);
xticklabels(ax3,{'' ''}); yticklabels(ax3,{'' '' ''})
ax3.LineWidth = linWidAx;
%
grid off;
%% Annotations
ax0 = x0+xLgd-0.101;
ay0 = 0.915+yLgd-0.57;
ady=0.054;
annotation('textbox',...
    [ax0+0.0005 0.900 0.02 0.02],...
    'String',{'*'},...
    'LineStyle','none',...
    'FontSize',fontSz*2,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');
annotation('textbox',...
    [ax0+0.003 0.8243 0.02 0.02],...
    'Color',[0 0 1],...
    'String','x',...
    'LineStyle','none',...
    'FontSize',fontSz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');
annotation('textbox',...
    [ax0+0.003 0.7489 0.02 0.02],...
    'Color',[0 0.5 0],...
    'String','x',...
    'LineStyle','none',...
    'FontSize',fontSz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');
annotation('textbox',...
    [ax0+0.003 0.6738 0.02 0.02],...
    'Color',[1 0 0],...
    'String','x',...
    'LineStyle','none',...
    'FontSize',fontSz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');
annotation('textbox',...
    [ax0-0.007 0.705 0.02 0.02],...
    'Color',[0 0.75 0.75],...
    'String','$\cdot$',...
    'LineStyle','none',...
    'FontSize',4*fontSz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off','interpreter','latex');
annotation('textbox',...
    [ax0-0.007 0.63 0.02 0.02],...
    'Color',[0.75 0 0.75],...
    'String','$\cdot$',...
    'LineStyle','none',...
    'FontSize',4*fontSz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off','interpreter','latex');
annotation('textbox',...
    [ax0-0.007 0.5575 0.02 0.02],...
    'Color',[0.75 0.75 0],...
    'String','$\cdot$',...
    'LineStyle','none',...
    'FontSize',4*fontSz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off','interpreter','latex');

%% Save
print(gcf,'../../results/sims3D.png','-dpng','-r500')
%print -depsc2 sims3D
