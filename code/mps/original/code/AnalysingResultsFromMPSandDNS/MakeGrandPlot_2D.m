%clear; clc;
restoredefaultpath;
clear restoredefaultpath;
addpath('../PDE_matfiles/General');
addpath('../plotTools/subplot_tight/subplot_tight');
addpath('../plotTools/colorbrewer');

load('extractedFields.mat'); load('processedFields.mat');

set(0,'DefaultLegendAutoUpdate','off')

% Placement
x0=0.13;
y0=0.81;
Y0=0.82;
ly=0.11;
lx=1.1*ly;
Lx=1.65*lx;
Ly=(340/390)*ly;
d0=0.128;
x1=0.79;
dy=0.77;
shift=0;%0.05;
dy0=1/7*dy+shift/2;
dy1=2/7*dy+shift/2;
dy2=3/7*dy+shift/2;
dy3=4/7*dy+shift;
dy4=5/7*dy+shift;
dy5=6/7*dy+shift;

fontsz = 10;
linWid=1.5;

labels1 =-1;

figure(1); clf
set(gcf, 'Units','points','Position',  [1920, 1080, 595, 430]) % A4 is 595 × 842 points. In arXiv template, default textwidth is 468 pt and text height is 648 pt.

%% Z-Vorticity plot

% DNS
subplot_tight(7,3,1,[0.00 0.00],'Position', [x0,y0,lx,ly]);
[ensMax, ensMin]=VisualiseFlowField_special(U250,V250 ,0,[1,1,1],fontsz,linWid,1,1);
ensMax=round(ensMax); ensMin=-ensMax;%round(ensMin);
caxis([ensMin ensMax]); camlight('left');
set(gca,'XAxisLocation','top','YAxisLocation','left');
text(labels1,0.5,'DNS', 'interpreter', 'latex','FontSize',fontsz)

subplot_tight(7,3,2,[0.00 0.00],'Position', [x0+d0*1,y0,lx,ly])
VisualiseFlowField_special(U750,V750,0,[1,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,3,[0.00 0.00],'Position', [x0+d0*2,y0,lx,ly])
VisualiseFlowField_special(U1250,V1250 ,0,[1,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,4,[0.00 0.00],'Position', [x0+d0*3,y0,lx,ly])
VisualiseFlowField_special(U1750,V1750 ,0,[1,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');
set(gca,'XAxisLocation','top','YAxisLocation','left');

% chi=118
subplot_tight(7,3,5,[0.00 0.00],'Position', [x0,y0-dy0,lx,ly])
ensMax=VisualiseFlowField_special(Uc118t250,Vc118t250 ,0,[0,1,0],fontsz,linWid,1,1);
ensMax=round(ensMax); ensMin=-ensMax;%round(ensMin);
caxis([ensMin ensMax]); camlight('left');
text(labels1,0.5,'$\chi{=}118$', 'interpreter', 'latex','FontSize',fontsz)
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,6,[0.00 0.00],'Position', [x0+d0*1,y0-dy0,lx,ly])
VisualiseFlowField_special(Uc118t750,Vc118t750,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,7,[0.00 0.00],'Position', [x0+d0*2,y0-dy0,lx,ly])
VisualiseFlowField_special(Uc118t1250,Vc118t1250,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,8,[0.00 0.00],'Position', [x0+d0*3,y0-dy0,lx,ly])
VisualiseFlowField_special(Uc118t1750,Vc118t1750,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

% chi=74
subplot_tight(7,3,9,[0.00 0.00],'Position', [x0,y0-dy1,lx,ly])
ensMax=VisualiseFlowField_special(Uc74t250,Vc74t250 ,0,[0,1,0],fontsz,linWid,1,1);
ensMax=round(ensMax); ensMin=-ensMax;%round(ensMin);
caxis([ensMin ensMax]); camlight('left');
text(labels1,0.5,'$\chi{=}74$', 'interpreter', 'latex','FontSize',fontsz)
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,10,[0.00 0.00],'Position', [x0+d0*1,y0-dy1,lx,ly])
VisualiseFlowField_special(Uc74t750,Vc74t750,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,11,[0.00 0.00],'Position', [x0+d0*2,y0-dy1,lx,ly])
VisualiseFlowField_special(Uc74t1250,Vc74t1250,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,12,[0.00 0.00],'Position', [x0+d0*3,y0-dy1,lx,ly])
VisualiseFlowField_special(Uc74t1750,Vc74t1750,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

% chi=33
subplot_tight(7,3,13,[0.00 0.00],'Position', [x0,y0-dy2,lx,ly])
ensMax=VisualiseFlowField_special(Uc33t250,Vc33t250 ,0,[0,1,0],fontsz,linWid,1,1);
ensMax=round(ensMax); ensMin=-ensMax;%round(ensMin);
caxis([ensMin ensMax]); camlight('left');
text(labels1,0.5,'$\chi{=}33$', 'interpreter', 'latex','FontSize',fontsz)
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,14,[0.00 0.00],'Position', [x0+d0*1,y0-dy2,lx,ly])
VisualiseFlowField_special(Uc33t750,Vc33t750,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,15,[0.00 0.00],'Position', [x0+d0*2,y0-dy2,lx,ly])
VisualiseFlowField_special(Uc33t1250,Vc33t1250,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,16,[0.00 0.00],'Position', [x0+d0*3,y0-dy2,lx,ly])
VisualiseFlowField_special(Uc33t1750,Vc33t1750,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

% 362^2
subplot_tight(7,3,17,[0.00 0.00],'Position', [x0,y0-dy3,lx,ly])
[ensMax, ensMin]=VisualiseFlowField_special(U250L8p5,V250L8p5 ,0,[0,1,0],fontsz,linWid,1,1);
ensMax=round(ensMax); ensMin=-ensMax;%round(ensMin);
caxis([ensMin ensMax]); camlight('left');
text(labels1,0.5,'$362^2\textnormal{ grid}$', 'interpreter', 'latex','FontSize',fontsz)
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,18,[0.00 0.00],'Position', [x0+d0*1,y0-dy3,lx,ly])
VisualiseFlowField_special(U750L8p5,V750L8p5,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,19,[0.00 0.00],'Position', [x0+d0*2,y0-dy3,lx,ly])
VisualiseFlowField_special(U1250L8p5,V1250L8p5 ,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,20,[0.00 0.00],'Position', [x0+d0*3,y0-dy3,lx,ly])
VisualiseFlowField_special(U1750L8p5,V1750L8p5 ,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

% 256^2
subplot_tight(7,3,21,[0.00 0.00],'Position', [x0,y0-dy4,lx,ly])
VisualiseFlowField_special(U250L8,V250L8,0,[0,1,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');
text(labels1,0.5,'$256^2\textnormal{ grid}$', 'interpreter', 'latex','FontSize',fontsz)
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,22,[0.00 0.00],'Position', [x0+d0*1,y0-dy4,lx,ly])
VisualiseFlowField_special(U750L8,V750L8,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,23,[0.00 0.00],'Position', [x0+d0*2,y0-dy4,lx,ly])
VisualiseFlowField_special(U1250L8,V1250L8,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,24,[0.00 0.00],'Position', [x0+d0*3,y0-dy4,lx,ly])
VisualiseFlowField_special(U1750L8,V1750L8,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

% 128^2
subplot_tight(7,3,25,[0.00 0.00],'Position', [x0,y0-dy5,lx,ly])
VisualiseFlowField_special(U250L7,V250L7,0,[0,1,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');
text(labels1,0.5,'$128^2\textnormal{ grid}$', 'interpreter', 'latex','FontSize',fontsz)
set(gca,'XAxisLocation','top','YAxisLocation','left');

subplot_tight(7,3,26,[0.00 0.00],'Position', [x0+d0*1,y0-dy5,lx,ly])
VisualiseFlowField_special(U750L7,V750L7,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

subplot_tight(7,3,27,[0.00 0.00],'Position', [x0+d0*2,y0-dy5,lx,ly])
VisualiseFlowField_special(U1250L7,V1250L7,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

ax1=subplot_tight(7,3,28,[0.00 0.00],'Position', [x0+d0*3,y0-dy5,lx,ly]);
VisualiseFlowField_special(U1750L7,V1750L7,0,[0,0,0],fontsz,linWid,1,1);
caxis([ensMin ensMax]); camlight('left');

%% Colorscale and colorbar for plot A
cmap=colormap(jet(256));

cmap( (128-round(5/97*256)):(128+round(5/97*256)),:)=1;
colormap(ax1,cmap);

%caxis([log10(10^-17) log10(1)]);
colorticks = [-97, -50, 0, 50, +97];
%ticklabels = {'$-100$','$-60$','$-20$','$0$','$20$','$60$','$100$'};
h=colorbar('Location','southoutside','Ticks',linspace(-97,97,11),'FontSize',fontsz,'Ticks',colorticks);%,'TickLabels',ticklabels);

set(h, 'Position', [x0 .08 .505 .05])
h.LineWidth=linWid;
h.TickLabelInterpreter = 'latex';
h.Label.Interpreter = 'latex';
h.Label.String = '$\nabla {\times} \mathit{\mathbf{V}}/T_0$';

%% Reynolds Stress 12
t=linspace(0,Tf,round(Tf/Ts+1));
N=length(insReStress12_MAT(:,1,1,1,1));
coord0=.2;
coord1=.803;

ResMax=max([insReStress12_MAT(:);insReStress12_TNT(:);]);
ResMin=min([insReStress12_MAT(:);insReStress12_TNT(:);]);

subplot_tight(7,1,1,[0.00 0.00],'Position',[x1,Y0,Lx,Ly])
pcolor(t,y(round(N*coord0):round(N*coord1)),(squeeze(insReStress12_MAT(round(N*coord0):round(N*coord1),4,1,:)))); shading interp;
ax=gca;
ax.LineWidth=linWid;set(gca,'layer','top');
ylabel('$y/L_{\mathrm{box}}$','interpreter','latex' ,'FontSize',fontsz);
xlabel('$t/T_0$','interpreter','latex' ,'FontSize',fontsz);
ax.FontSize=fontsz;
%xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );
grid on;
c0=caxis;%camlight
caxis([ResMin ResMax]);
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
text(-1.15,0.5,'DNS', 'interpreter', 'latex','FontSize',fontsz)

tau12=(squeeze(insReStress12_MAT(round(N*coord0):round(N*coord1),4,1,:)));
tau12_max=max(tau12(:));
tau12_min=min(tau12(:));

cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);

subplot_tight(7,1,2,[0.00 0.00],'Position',[x1,Y0-dy0,Lx,Ly])
pcolor(t,y(round(N*coord0):round(N*coord1)),(squeeze(insReStress12_TNT(round(N*coord0):round(N*coord1),4,1,3,:)))); shading interp;
c2= caxis;%camlight
caxis([ResMin ResMax])
ax=gca;
ax.LineWidth=linWid; set(gca,'layer','top');%set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );yticks([0.2 .5 .8] )
grid on;
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'Xticklabel',[]);ax.FontSize=fontsz;
text(-1.15,0.5,'$\chi{=}118$', 'interpreter', 'latex','FontSize',fontsz)
 
cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);

subplot_tight(7,1,3,[0.00 0.00],'Position',[x1,Y0-dy1,Lx,Ly])
pcolor(t,y(round(N*coord0):round(N*coord1)),(squeeze(insReStress12_TNT(round(N*coord0):round(N*coord1),4,1,2,:)))); shading interp;
%text(0.7,1.1,'$\chi=60$ ($\sim$1:16)','interpreter','latex' , 'FontSize',30);
c4= caxis;%camlight
caxis([ResMin ResMax])
ax=gca;
ax.LineWidth=linWid; set(gca,'layer','top');%set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );yticks([0.2 .5 .8] )
 grid on;
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'Xticklabel',[]);ax.FontSize=fontsz;
text(-1.15,0.5,'$\chi{=}74$', 'interpreter', 'latex','FontSize',fontsz)

cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);
 
subplot_tight(7,1,4,[0.00 0.00],'Position',[x1,Y0-dy2,Lx,Ly])
pcolor(t,y(round(N*coord0):round(N*coord1)),(squeeze(insReStress12_TNT(round(N*coord0):round(N*coord1),4,1,1,:)))); shading interp;
c6= caxis;%camlight
caxis([ResMin ResMax])
ax=gca;
ax.LineWidth=linWid; set(gca,'layer','top');
%set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );yticks([0.2 .5 .8] )
 grid on;
ax=gca;
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'Xticklabel',[]);ax.FontSize=fontsz;
text(-1.15,0.5,'$\chi{=}33$', 'interpreter', 'latex','FontSize',fontsz)

cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);

subplot_tight(7,1,5,[0.00 0.00],'Position',[x1,Y0-dy3,Lx,Ly])
pcolor(t,linspace(y(round(N*coord0)),y(round(N*coord1)),length(round(362*coord0):round(362*coord1))),(squeeze(insReStress12_MAT(round(362*coord0):round(362*coord1),3,1,:)))); shading interp;
c1= caxis;%camlight
caxis([ResMin ResMax]); grid on;
ax=gca;
ax.LineWidth=linWid; set(gca,'layer','top');%set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );yticks([0.2 .5 .8] )
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'Xticklabel',[]);ax.FontSize=fontsz;
text(-1.15,0.5,'$362^2 \textnormal{ grid}$', 'interpreter', 'latex','FontSize',fontsz)

cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);

subplot_tight(7,1,6,[0.00 0.00],'Position',[x1,Y0-dy4,Lx,Ly])
pcolor(t,linspace(y(round(N*coord0)),y(round(N*coord1)),length(round(256*coord0):round(256*coord1))),(squeeze(insReStress12_MAT(round(256*coord0):round(256*coord1),2,1,:)))); shading interp;
c3= caxis;%camlight
caxis([ResMin ResMax])
ax=gca;
ax.LineWidth=linWid; set(gca,'layer','top');%set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );yticks([0.2 .5 .8] )
grid on;
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'Xticklabel',[]);ax.FontSize=fontsz;
text(-1.15,0.5,'$256^2 \textnormal{ grid}$', 'interpreter', 'latex','FontSize',fontsz)
 
cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);

ax=subplot_tight(7,1,7,[0.00 0.00],'Position',[x1,Y0-dy5,Lx,Ly]);
pcolor(t,linspace(y(round(N*coord0)),y(round(N*coord1)),length(round(128*coord0):round(128*coord1))),(squeeze(insReStress12_MAT(round(128*coord0):round(128*coord1),1,1,:)))); shading interp;
c5= caxis;
%camlight('left');
caxis([ResMin ResMax])
ax.LineWidth=linWid; set(gca,'layer','top');%set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);xticks([0 0.25 .5 .75 1 1.25 1.5 1.75 2] );yticks([0.2 .5 .8] )
grid on;
xticks([0 .5 1 1.5 2] );
yticks([0.3 .5 .7] ); 
set(gca,'XAxisLocation','top','YAxisLocation','left');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'Xticklabel',[]);ax.FontSize=fontsz;
text(-1.15,0.5,'$128^2 \textnormal{ grid}$', 'interpreter', 'latex','FontSize',fontsz)

cmap=colormap(jet(256));
cmap( (128-round(-0.0025/tau12_min*256)):(128+round(0.001/tau12_max*256)),:)=1;
colormap(ax,cmap);

%% Colorscale and colorbar for plot B

% %cFull=[min([c0 c1 c2 c3 c4 c5 c6]), max([c0 c1 c2 c3 c4 c5 c6])];
% 
% NumColors= 256;
% wbandL=5;
% wbandR=4;
% 
% cmap=colormap(jet(NumColors));
% cmap(round(NumColors/2-NumColors/39-1*wbandL):round(NumColors/2-NumColors/39+wbandR),:)=1*ones(size(cmap(round(NumColors/2-NumColors/35-1*wbandL):round(NumColors/2-NumColors/35+wbandR),:)));
set(gcf,'color','white')
%
% colormap(cmap);

%caxis([log10(10^-17) log10(1)]);
%colorticks = [-97, -50, 0, 50, +97];
ticklabels = {'$-100$','$-60$','$-20$','$0$','$20$','$60$','$100$'};
h=colorbar('Location','southoutside');%,'Ticks',linspace(-1,1,11),'FontSize',fontsz);%,'Ticks',colorticks);%,'TickLabels',ticklabels);

set(h, 'Position', [x1 .08 Lx .05])
h.LineWidth=linWid;
h.TickLabelInterpreter = 'latex';
h.Label.Interpreter = 'latex';
h.Label.String = '$\tau_{12}/(L_{\textrm{box}}/T_0)^2$';

%% Annonation
ax0 = 0.66;


annotation('line',[ax0 ax0],...
    [1 0],'LineWidth',1.5);

annotation('textbox',...
    [0.008 0.98 0.02 0.02],...
    'String','\textbf{a}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation('textbox',...
    [ax0+0.008 0.98 0.02 0.02],...
    'String','\textbf{b}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',1.2*fontsz,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

%% Save
%set(gcf, 'Renderer', 'Painters');
%print -depsc2 Sims2D

print(gcf,'../../results/sims2D.png','-dpng','-r500')