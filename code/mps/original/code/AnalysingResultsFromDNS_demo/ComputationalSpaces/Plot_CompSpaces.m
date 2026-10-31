% Plot M, W and S in 2D & 3D.
% In 2D: W is set on 128^2 grid, M is set at chi_max. 4 DNS snapshots are plotted 
% @ t=0.25,0.75,1.25,1.75 at a certain accuracy.
% In 3D: W is set on 64^3 grid, M is set at chi_max. 4 DNS snapshots are plotted 
% @ t=0.2,0.8,1.4,2 at a certain accuracy.

%% First addpath
clear all;
addpath('../../plotTools/subplot_tight/subplot_tight');
addpath('../../plotTools/legendflex');

%% Plot parameters
set(gcf,'color','w');

% Placement
xp0= 0.07;
yp0= 0.088;
lx = 0.3;
ly = 0.34;
%
xp1 = 0.47;
yp1= 0.59;
lx1 = 0.33;
ly1=0.35;
%
xp2 = xp1;
yp2=yp0;
lx2 = lx1;
ly2=0.35;
%
linWid = 1.5;

linWid_ss=1;
mrkSz=4;
fontsz=10;

count = 3;

% Legends
%chi vs re
legends_chiRe = {'\textnormal{TGV} (3-D)','$\chi_{99}{\sim}\textnormal{Re}^{0.71}$','\textnormal{TDJ} (2-D)'};
lgdChiRe_pos= [-95 -250];
lgdChiRe_width = [0.0,0.0,.1,.1];

%2D comp space
legends2D = {'$\textnormal{TDJ},t/T_0{=}0.25 \textnormal{    }$',...
    '$\textnormal{TDJ},t/T_0{=}0.75$','$\textnormal{TDJ},t/T_0{=}1.25$','$\textnormal{TDJ},t/T_0{=}1.75$','$\mathcal{D}, 256^2 \textnormal{ grid}$','$\mathcal{W}, 256^2 \textnormal{ grid}$','$\mathcal{M},\chi_{99}{=}20$'};
lgd2D_pos = [+125 +0];
lgd2D_width = [0.0,0.0,.3,.3];

%3D comp space
legends3D = {'$\textnormal{TGV},t/T_0{=}0.2 \textnormal{    }$',...
    '$\textnormal{TGV},t/T_0{=}0.8$','$\textnormal{TGV},t/T_0{=}1.4$','$\textnormal{TGV},t/T_0{=}2$','$\mathcal{D}, 64^3 \textnormal{ grid}$','$\mathcal{W}, 64^3 \textnormal{ grid}$','$\mathcal{M},\chi_{99}{=}41$'};%,'$\mathcal{W}, 16^3 \textnormal{ grid}$'};
lgd3D_pos = [+125 0];
lgd3D_width = [0.0,0.0,.3,.3];

%% First plot the chi vs Re connection
clf; 
load('IncompNS_ChisRe');


gcf=figure(1);
set(gcf, 'Units','points','Position',  [0, 0, 468, 324]) % A4 is 595 × 842 points or 210 × 297 millimeters or 8.27 × 11.69 inches. In arXiv template, default textwidth is 468 pt and text height is 648 pt.


subplot_tight(2,2,2,[0.00 0.00],'Position',[xp0,yp0,lx,ly])

currPlot(1)=plot(Re3D,chis3D_log2,'d','MarkerSize',mrkSz,'LineWidth',linWid,'color', [0,0,1]);hold on;
currPlot(2)=plot(Re3D,nikEstimate3D_log2+(chis3D_log2(2)-nikEstimate3D_log2(2)),'-','LineWidth',linWid,'color', [0,0,1]);
currPlot(3)=plot(Re2D,chis2D_log2,'.','MarkerSize',4*mrkSz,'LineWidth',linWid,'color', [0,0.5,0]);

%Plot design
yticks(0:2:8);
ylim([0 8]);
xticks([10,100,1000]);

axesH = gca;
set(axesH, 'YMinorTick','on', 'YMinorGrid','on')
axesH.YAxis.MinorTickValues = [1,3,5,7];
set(gca, 'XScale', 'log')

[lgd,objh] = legendflex(currPlot,legends_chiRe, 'xscale',0.2, 'ref', gcf, 'anchor', [2,2], 'buffer',lgdChiRe_pos, 'nrow',3,'ncol',1, 'FontSize',fontsz,'Position', lgdChiRe_width, 'Interpreter','latex','box','off');
lgd.FontSize = fontsz;
lineh = findobj(objh,'type','line');
set(lineh,'LineWidth',linWid)
lgd.LineWidth = linWid;

yticklabels({'$1$','$4$','$4^2$','$4^3$','$4^4$'});
xlabel('$\textnormal{Re}$','interpreter','latex');
ylabel('$\chi_{99}$','interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = fontsz;
ax.LineWidth = linWid; % Make tick marks thicker.
grid off;

%% Now plot 2D
load('IncompNS2D/MWS_e2to12_out');%fprintf('\n2D Accuracy is: %.2s',desiredAcc(count));
%
h=subplot_tight(2,2,3,[0.00 0.00],'Position',[xp1,yp1,lx1,ly1]);
hold on;
%
chi_max = 20;
x0=2.35; x1=5.65;
%l_W=log2(64);
%
%First plot S
AS=area(1:L-1,log2(chiMax));
set(AS,'FaceColor',1*[1 1 1]);
AS.FaceAlpha = 0;

AS.LineWidth = linWid; 
AS.LineStyle = '--';
set(AS,'EdgeColor',[0,0,0]);
set(AS,'EdgeAlpha',.8);

%Now plot M at appropriate chi_max
Y_mps = [2,4,log2(chi_max),log2(chi_max),4,2];
X_mps = [1,2,x0,x1,6,7];
AMPS=area(X_mps,Y_mps);
set(AMPS,'FaceColor',1*[0 0.4470 0.7410]);
AMPS.FaceAlpha = .2;
%
AMPS.LineWidth = linWid; 
AMPS.LineStyle = ':';
set(AMPS,'EdgeColor',1*[0 0.4470 0.7410]);
set(AMPS,'EdgeAlpha',0.5);

%Now plot W
l_W=8;
chiMax_8 = [4,4^2,4^3,4^4,4^3,4^2,4^1];
AW=area(1:(l_W-1),log2(chiMax_8));
set(AW,'FaceColor',0*[1,1,1]);
AW.FaceAlpha = .2;
%
AW.LineWidth = linWid; 
AW.LineStyle = '-.';
set(AW,'EdgeColor',0*[1,1,1]);
set(AW,'EdgeAlpha',0.5);

%And finally plot DNS flow snapshots
%t=.25
%plot(1:L-1,log2(chi250),'linewidth',linWid,'color',[0.8500 0.3250 0.0980],'LineStyle','-');
A1=area(1:L-1,log2(chi250(:,count)));
set(A1,'FaceColor',[0.8500 0.3250 0.0980]);
A1.FaceAlpha =.3;

A1.LineWidth = linWid_ss; 
A1.LineStyle = '-';
set(A1,'EdgeColor',[0,0,1]);
set(A1,'EdgeAlpha',.5);

%t=.75
%plot(1:L-1,log2(chi750),'linewidth',linWid,'color',[0.8500 0.3250 0.0980],'LineStyle','-');
A2=area(1:L-1,log2(chi750(:,count)));
set(A2,'FaceColor',[0.8500 0.3250 0.0980]); A2.EdgeAlpha = 0;
A2.FaceAlpha = .3;

A2.LineWidth = linWid_ss; 
A2.LineStyle = '-';
set(A2,'EdgeColor',[0,0.5,0]);
set(A2,'EdgeAlpha',.5);

%t=1.25
%plot(1:L-1,log2(chi1250),'linewidth',linWid,'color',[0.8500 0.3250 0.0980],'LineStyle','-');
A3=area(1:L-1,log2(chi1250(:,count)));
set(A3,'FaceColor',[0.8500 0.3250 0.0980]);
A3.FaceAlpha = .3;

A3.LineWidth = linWid_ss; 
A3.LineStyle = '-';
set(A3,'EdgeColor',[1,0,0]);
set(A3,'EdgeAlpha',.5);

%t=1.75
%plot(1:L-1,log2(chi1750),'linewidth',linWid,'color',[0.8500 0.3250 0.0980],'LineStyle','-');
A4=area(1:L-1,log2(chi1750(:,count)));
set(A4,'FaceColor',1*[0.8500 0.3250 0.0980]);
A4.FaceAlpha = .3;

A4.LineWidth = linWid_ss; 
A4.LineStyle = '-';
set(A4,'EdgeColor',[0,0.75,.75]);
set(A4,'EdgeAlpha',.5);
%
[lgd,h_legend] = legendflex([A1,A2,A3,A4,AS,AW,AMPS],legends2D, 'xscale',0.4,'ref', h, 'anchor', [2,2], 'buffer',lgd2D_pos, 'nrow',7,'ncol',1, 'FontSize',fontsz,'Position', lgd2D_width,'interpreter','latex','Box','off');
lgd.LineWidth = linWid;
PatchInLegend = findobj(h_legend, 'type', 'patch');
set(PatchInLegend(1), 'FaceAlpha', 0.2);set(PatchInLegend(1), 'EdgeAlpha', 0.5);
set(PatchInLegend(2), 'FaceAlpha', 0.3);set(PatchInLegend(2), 'EdgeAlpha', 0.5);
set(PatchInLegend(3), 'FaceAlpha', 0.3);set(PatchInLegend(3), 'EdgeAlpha', 0.5);
set(PatchInLegend(4), 'FaceAlpha', 0.3);set(PatchInLegend(4), 'EdgeAlpha', 0.5);
set(PatchInLegend(5), 'FaceAlpha', 0.3);set(PatchInLegend(5), 'EdgeAlpha', 0.5);
set(PatchInLegend(6), 'FaceAlpha', 0.3);set(PatchInLegend(6), 'EdgeAlpha', 0.5);
set(PatchInLegend(7), 'FaceAlpha', 0.3);set(PatchInLegend(7), 'EdgeAlpha', 0.5);

%
% 2D Plot settings
xlim([1 9]);
xticks(1:9);
%xticklabels({'$2^9$','$2^8$','$2^7$','$2^6$','$2^5$','$2^4$','$2^3$','$2^2$','$2$'});
xticklabels({'$1$','$2$','$3$','$4$','$5$','$6$','$7$','$8$','$9$'});
ylim([0 10]);
yticks(0:2:10);
yticklabels({'$1$','$4$','$4^2$','$4^3$','$4^4$','$4^5$'});
%
axesH = gca;
set(axesH, 'YMinorTick','on', 'YMinorGrid','on')
axesH.YAxis.MinorTickValues = 1:2:log2(N)-1;
%
xlabel('$n$-th Bipartition','interpreter','latex');
%xlabel('$\frac{\Delta x_{\S^{\mathrm{max}}_{99}{rough}}}{\Delta x_{\S^{\mathrm{max}}_{99}{fine}}}$','interpreter','latex');
ylabel('$d(n)$','interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = fontsz;
ax.LineWidth = linWid; % Make tick marks thicker.
grid off;

%% Then plot 3D

load('IncompNS3D/MWS_e2to12_out');%fprintf('\n3D Accuracy is: %.2s\n\n',desiredAcc(count));
%
h=subplot_tight(2,2,4,[0.00 0.00],'Position',[xp2,yp2,lx2,ly2]);
hold on;
%
chi_max = 41;
x0=1.66; x1=4.32;
l_W=log2(32);
%First plot S
AS=area(1:L-1,log2(chiMax));
set(AS,'FaceColor',1*[1 1 1]);
AS.FaceAlpha = 0;

AS.LineWidth = linWid; 
AS.LineStyle = '--';
set(AS,'EdgeColor',[0,0,0]);
set(AS,'EdgeAlpha',.8);

%Now plot M at appropriate chi_max
Y_mps = [3,log2(chi_max),log2(chi_max),3];
X_mps = [1,x0,x1,5];
AMPS=area(X_mps,Y_mps);
set(AMPS,'FaceColor',1*[0 0.4470 0.7410]);
AMPS.FaceAlpha = .2;
%
AMPS.LineWidth = linWid; 
AMPS.LineStyle = ':';
set(AMPS,'EdgeColor',1*[0 0.4470 0.7410]);
set(AMPS,'EdgeAlpha',0.5);

%Now plot W
l_W=6;
chiMax_6 = [8,8^2,8^3,8^2,8];
AW=area(1:(l_W-1),log2(chiMax_6));
set(AW,'FaceColor',0*[1,1,1]); AW.EdgeAlpha = 0;
AW.FaceAlpha = .2;
%
AW.LineWidth = linWid; 
AW.LineStyle = '-.';
set(AW,'EdgeColor',0*[1,1,1]);
set(AW,'EdgeAlpha',0.5);

%And finally plot DNS flow snapshots
%t=.25
A1=area(1:L-1,log2(chi200(:,count)));
set(A1,'FaceColor',[0.8500 0.3250 0.0980]);
A1.FaceAlpha =.3;

A1.LineWidth = linWid_ss; 
A1.LineStyle = '-';
set(A1,'EdgeColor',[0,0,1]);
set(A1,'EdgeAlpha',.5);

%t=.75
A2=area(1:L-1,log2(chi800(:,count)));
set(A2,'FaceColor',[0.8500 0.3250 0.0980]); A2.EdgeAlpha = 0;
A2.FaceAlpha = .3;

A2.LineWidth = linWid_ss; 
A2.LineStyle = '-';
set(A2,'EdgeColor',[0,0.5,0]);
set(A2,'EdgeAlpha',.5);

%t=1.25
A3=area(1:L-1,log2(chi1400(:,count)));
set(A3,'FaceColor',[0.8500 0.3250 0.0980]);
A3.FaceAlpha = .3;

A3.LineWidth = linWid_ss; 
A3.LineStyle = '-';
set(A3,'EdgeColor',[1,0,0]);
set(A3,'EdgeAlpha',.5);

%t=1.75
%plot(1:L-1,log2(chi1750),'linewidth',linWid,'color',[0.8500 0.3250 0.0980],'LineStyle','-');
A4=area(1:L-1,log2(chi2000(:,count)));
set(A4,'FaceColor',1*[0.8500 0.3250 0.0980]);
A4.FaceAlpha = .3;

A4.LineWidth = linWid_ss; 
A4.LineStyle = '-';
set(A4,'EdgeColor',[0,0.75,.75]);
set(A4,'EdgeAlpha',.5);
%
%
[lgd,h_legend] = legendflex([A1,A2,A3,A4,AS,AW,AMPS],legends3D, 'xscale',0.4, 'ref', h, 'anchor', [2,2], 'buffer',lgd3D_pos, 'nrow',7,'ncol',1, 'FontSize',fontsz,'Position', lgd3D_width,'interpreter','latex','Box','off');
PatchInLegend = findobj(h_legend, 'type', 'patch');
set(PatchInLegend(1), 'FaceAlpha', 0.2);set(PatchInLegend(1), 'EdgeAlpha', 0.5);
set(PatchInLegend(2), 'FaceAlpha', 0.3);set(PatchInLegend(2), 'EdgeAlpha', 0.5);
set(PatchInLegend(3), 'FaceAlpha', 0.3);set(PatchInLegend(3), 'EdgeAlpha', 0.5);
set(PatchInLegend(4), 'FaceAlpha', 0.3);set(PatchInLegend(4), 'EdgeAlpha', 0.5);
set(PatchInLegend(5), 'FaceAlpha', 0.3);set(PatchInLegend(5), 'EdgeAlpha', 0.5);
set(PatchInLegend(6), 'FaceAlpha', 0.3);set(PatchInLegend(6), 'EdgeAlpha', 0.5);
set(PatchInLegend(7), 'FaceAlpha', 0.3);set(PatchInLegend(7), 'EdgeAlpha', 0.5);

% 3D Plot settings
xticks(1:7);
xlim([1 7]);
%xticklabels({'$2^7$','$2^6$','$2^5$','$2^4$','$2^3$','$2^2$','$2$'});
xticklabels({'$1$','$2$','$3$','$4$','$5$','$6$','$7$'});
yticks(0:3:12);
ylim([0 12]);
%
axesH = gca;
set(axesH, 'YMinorTick','on', 'YMinorGrid','on')
axesH.YAxis.MinorTickValues = [1,2,4,5,7,8,10,11];
%
%set(gca, 'YMinorTick','on', 'YMinorGrid','on')
yticklabels({'$1$','$8$','$8^2$','$8^3$','$8^4$'});
%xlabel('$\frac{\Delta x_{\S^{\mathrm{max}}_{99}{rough}}}{\Delta x_{\S^{\mathrm{max}}_{99}{fine}}}$','interpreter','latex');
xlabel('$n$-th Bipartition','interpreter','latex');
ylabel('$d(n)$','interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = fontsz;
ax.LineWidth = linWid; % Make tick marks thicker.
grid off;

%% Annotations
annotation('textbox',...
    [0 .995 0 0.02],...
    'String','\textbf{a}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation('textbox',...
    [0.45 0.995 0.02 0.02],...
    'String','\textbf{b}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation('textbox',...
    [0 0.5 0 0.02],...
    'String','\textbf{d}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation('textbox',...
    [0.45 0.5 0.02 0.02],...
    'String','\textbf{c}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

%annotation('line',[0.418 0.418], [1 0.51],'LineWidth',6);
%annotation('line',[0 1], [0.51 0.51],'LineWidth',3);
%annotation('line',[0.59 .59], [0.513 0],'LineWidth',3);
% 
%print -depsc2 CompSpaces
print(gcf,'../../../results/CompSpaces_demo.png','-dpng','-r500')