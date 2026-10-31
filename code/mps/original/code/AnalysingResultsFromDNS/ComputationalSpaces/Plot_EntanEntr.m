%% First clear and add paths
clear;
addpath('.../../plotTools/subplot_tight/subplot_tight');
addpath('../../plotTools/legendflex');

%% Plot parameters
set(gcf,'color','w');

% Placement
x0= 0.08;
y0= 0.146;

lx = 0.4;
ly = 0.8;

x1 = x0+lx + 0.1;
lx1 = lx;

figure(1); clf
linWid = 1.5;
mrkSz=4;
fontsz=10;

set(gcf, 'Units','points','Position',  [1920, 1080, 468, 324]) % A4 is 595 × 842 points. In arXiv template, default textwidth is 468 pt and text height is 648 pt.

%
lgd_pos = [-58 -21];
lgd_pos3D = [+172 -21];
lgd_width = [0,0,.01,.1];


%% Entanglement Entropy 2D
%
load('IncompNS2D/MWS_e2to12_out');
%
legends = {'$u_1,t/T_0{=}0.25\;$','$u_1 , t/T_0{=}0.75$','$u_1,t/T_0{=}1.25$','$u_1 ,t/T_0{=}1.75$','$u_2,t/T_0{=}0.25$','$u_2 , t/T_0{=}0.75$','$u_2,t/T_0{=}1.25$','$u_2 ,t/T_0{=}1.75$'};
subplot_tight(1,2,1,[0.00 0.00],'Position',[x0,y0,lx,ly])

biparts = 1:9;

%U
currPlot(1)=plot(biparts,entanEntrU250(1,:),'LineStyle','-','LineWidth',linWid,'color', [0,0,1]); hold on;
currPlot(2)=plot(biparts,entanEntrU750(1,:),'LineStyle','--','LineWidth',linWid,'color', [0,0,1]);
currPlot(3)=plot(biparts,entanEntrU1250(1,:),'LineStyle','-.','LineWidth',linWid,'color', [0,0,1]);
currPlot(4)=plot(biparts,entanEntrU1750(1,:),'LineStyle',':','LineWidth',linWid,'color', [0,0,1]);
%V
currPlot(5)=plot(biparts,entanEntrV250(1,:),'LineStyle','-','LineWidth',linWid,'color', [0,0.5,0]); hold on;
currPlot(6)=plot(biparts,entanEntrV750(1,:),'LineStyle','--','LineWidth',linWid,'color', [0,0.5,0]);
currPlot(7)=plot(biparts,entanEntrV1250(1,:),'LineStyle','-.','LineWidth',linWid,'color', [0,0.5,0]);
currPlot(8)=plot(biparts,entanEntrV1750(1,:),'LineStyle',':','LineWidth',linWid,'color', [0,0.5,0]);

%Axis and labels
maxEntanEntrInPractice = max([entanEntrU250(1,:),entanEntrU750(1,:),entanEntrU1250(1,:),entanEntrU1750(1,:),...
    entanEntrV250(1,:),entanEntrV750(1,:),entanEntrV1250(1,:),entanEntrV1750(1,:)]);
%
xticks(biparts);
yticks(linspace(0,ceil(maxEntanEntrInPractice),7));
xlim([1 max(biparts)]);
ylim([0,ceil(maxEntanEntrInPractice)]);
%
%xticklabels({'$2^9$','$2^8$','$2^7$','$2^6$','$2^5$','$2^4$','$2^3$','$2^2$','$2$'});
xticklabels({'$1$','$2$','$3$','$4$','$5$','$6$','$7$','$8$','$9$'});
%yticklabels({'$$'});
set(gca,'TickLabelInterpreter', 'latex');
ax=gca;
ax.LineWidth = linWid; set(gca,'layer','top');
ax.FontSize = fontsz;
%ylabel('$H[\Delta x_{\mathrm{rough}}/\Delta x_{\mathrm{fine}}]$','interpreter','latex','fontsize',fontsz);
%xlabel('$\Delta x_{\mathrm{rough}}/\Delta x_{\mathrm{fine}}$','interpreter','latex');
ylabel('$H(n,t)$','interpreter','latex','fontsize',fontsz);
xlabel('$n$-th Bipartition','interpreter','latex');

%[lgd,objh] = legend((currPlot),legends,'Interpreter','latex','Position',[x0+0.322 y0+.78*ly 0 0]);
[lgd,objh] = legendflex(currPlot,legends, 'xscale',0.7,'ref', gcf, 'anchor', [2,2], 'buffer',lgd_pos, 'nrow',8,'ncol',1, 'FontSize',fontsz,'Position', lgd_width, 'Interpreter','latex');
lineh = findobj(objh,'type','line');
set(lineh,'LineWidth',linWid)
lgd.LineWidth = linWid;
grid off;
%

%% Entanglement Entropy 3D
%
load('IncompNS3D/MWS_e2to12_out');
%
legends = {'$u_1,u_2; t/T_0{=}0.2\;$','$u_1,u_2; t/T_0{=}0.8$','$u_1,u_2;t/T_0{=}1.4$','$u_1,u_2;t/T_0{=}2$','$u_3,t/T_0{=}0.2$','$u_3,t/T_0{=}0.8$','$u_3,t/T_0{=}1.4$','$u_3,t/T_0{=}2$'};

subplot_tight(1,2,2,[0.00 0.00],'Position',[x1,y0,lx1,ly])

biparts = 1:7;

% U, V
currPlot(1)=plot(biparts,entanEntrU200(1,:),'LineStyle','-','LineWidth',linWid,'color', [0,0,1]); hold on;
currPlot(2)=plot(biparts,entanEntrU800(1,:),'LineStyle','--','LineWidth',linWid,'color', [0,0,1]);
currPlot(3)=plot(biparts,entanEntrU1400(1,:),'LineStyle','-.','LineWidth',linWid,'color', [0,0,1]);
currPlot(4)=plot(biparts,entanEntrU2000(1,:),'LineStyle',':','LineWidth',linWid,'color', [0,0,1]);
%W
currPlot(5)=plot(biparts,entanEntrW200(1,:),'LineStyle','-','LineWidth',linWid,'color', [0,0.5,0]); hold on;
currPlot(6)=plot(biparts,entanEntrW800(1,:),'LineStyle','--','LineWidth',linWid,'color', [0,0.5,0]);
currPlot(7)=plot(biparts,entanEntrW1400(1,:),'LineStyle','-.','LineWidth',linWid,'color', [0,0.5,0]);
currPlot(8)=plot(biparts,entanEntrW2000(1,:),'LineStyle',':','LineWidth',linWid,'color', [0,0.5,0]);

%Axis and labels
maxEntanEntrInPractice = max([entanEntrU200(1,:),entanEntrU800(1,:),entanEntrU1400(1,:),entanEntrU2000(1,:),...
    entanEntrV200(1,:),entanEntrV800(1,:),entanEntrV1400(1,:),entanEntrV2000(1,:),entanEntrW200(1,:),entanEntrW800(1,:),...
    entanEntrW1400(1,:),entanEntrW2000(1,:)]);
%
xticks(biparts);
yticks(linspace(0,ceil(maxEntanEntrInPractice),6));
xlim([1 max(biparts)]);
ylim([0,ceil(maxEntanEntrInPractice)]);
%
%xticklabels({'$2^7$','$2^6$','$2^5$','$2^4$','$2^3$','$2^2$','$2$'});
xticklabels({'$1$','$2$','$3$','$4$','$5$','$6$','$7$'});
%yticklabels({'$$'});
set(gca,'TickLabelInterpreter', 'latex');
ax=gca;
ax.LineWidth = linWid; set(gca,'layer','top');
ax.FontSize = fontsz;
%ylabel('$H[\Delta x_{\mathrm{rough}}/\Delta x_{\mathrm{fine}}]$','interpreter','latex','fontsize',fontsz);
%xlabel('$\Delta x_{\mathrm{rough}}/\Delta x_{\mathrm{fine}}$','interpreter','latex');
ylabel('$H(n,t)$','interpreter','latex','fontsize',fontsz);
xlabel('$n$-th Bipartition','interpreter','latex');
[lgd,objh] = legendflex(currPlot,legends, 'xscale',0.7, 'ref', gcf, 'anchor', [2,2], 'buffer',lgd_pos3D, 'nrow',8,'ncol',1, 'FontSize',fontsz,'Position', lgd_width, 'Interpreter','latex');
lgd.FontSize = fontsz;
lineh = findobj(objh,'type','line');
set(lineh,'LineWidth',linWid)
lgd.LineWidth = linWid;
grid off;
%

%% Annotations
annotation('textbox',...
    [0.015 0.97 0.02 0.02],...
    'String','\textbf{a}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation('textbox',...
    [0.515 0.97 0.02 0.02],...
    'String','\textbf{b}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation('line',[0.51 0.51],...
    [1 0],'LineWidth',1.5);

%% Save
%print -depsc2 ../../results/EntanEntr
print(gcf,'../../../results/EntanEntr.png','-dpng','-r250')