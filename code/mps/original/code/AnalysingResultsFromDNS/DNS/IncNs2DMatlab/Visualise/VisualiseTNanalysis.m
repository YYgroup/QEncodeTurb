%Visualises various results from simulate2DFlow runs

addpath('/Volumes/ARBEID/OxfordFysikkJobb/NSProject/CompareTNTMATLAB');

plot_L = log2(size(U)); plot_L=plot_L(1);
plot_tvec = linspace(0,t,length(U_SingMatMax(1,:)));
D_max = min(length(U_errT(:,1)),1024);
plot_chivec = 0:(D_max-1);

%U_errT = U_errT(:,end);
%V_errT = V_errT(:,end);

%% Plot Entanglement Entropy as a function of time
figure(1)
clf;
subplot(1,2,1)
plot(plot_tvec,U_SingMatMax(1,:),'b','LineWidth',3); hold on;
plot(plot_tvec,V_SingMatMax(1,:),'r','LineWidth',3);
axis([0 t 0 max(max(U_SingMatMax(2,:),V_SingMatMax(2,:))) ])
title(sprintf('Maximum Entanglement Entropy of flow field\nRe = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
xlabel('t (seconds)');
ylabel('S_{max}(t)');
legend(' V_{x}',' V_{y}');


set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')

subplot(1,2,2)
plot(plot_tvec,U_SingMatMax(2,:),'b','LineWidth',3); hold on;
plot(plot_tvec,V_SingMatMax(2,:),'r','LineWidth',3);
axis([0 t 1 (plot_L-1) ])
title(sprintf('Bond at which Entanglement Entropy is Maximum of flow field\nRe = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
xlabel('t (seconds)');
ylabel('bond #');
legend(' V_{x}',' V_{y}');

set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')

%% Plot Error in Entanglement Entropy as a function of time
%Calculate fraction of energy uncaptured
norm = plot_L-1;
fracEnUncap = ( U_errT + V_errT )./( sum(U_errT,1) + sum(V_errT,1) );

log10compRatio = zeros(1,length(plot_chivec));
for i = 1:length(plot_chivec)
    log10compRatio(i) = log10(memComp(4,plot_chivec(i),plot_L));
end

%log10compRatio(1) = -16;

figure(2)
colormap('hsv')
%clf;
%subplot(1,2,1)
surf(plot_tvec,log10compRatio,log10( fracEnUncap ));colorbar; view(180,90);  shading interp; 
%plot(log10compRatio,log10(fracEnUncap(:,end)));

%axis([0 t 0 D_max 0 1 ])
title(sprintf('Energy uncaptured by MPS representation by compression ratio, across t.\nRe = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
xlabel('t (seconds)');
ylabel('Log10 Compression ratio (function of \chi)');
%zlabel('Log10 Error');

set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')

% subplot(1,2,2)
% norm = plot_L-1;
% surf(plot_tvec,plot_chivec,log10(V_errT(plot_chivec+1,:)./norm) ); view(120,20); shading interp;
% %axis([0 t 0 D_max 0 1 ])
% title(sprintf('Normalised Cumulative sum of Log-Singular Values squared of V_{y}\n up to each chi, across t. Re = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
% xlabel('t (seconds)');
% ylabel('\chi');
% zlabel('Normalised CumulativeSum(\chi,t)');
% 
% set(gca,'FontSize',24,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')