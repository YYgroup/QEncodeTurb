%Visualises various results from simulate2DFlow runs

plot_L = log2(size(U)); plot_L=plot_L(1);
plot_tvec = linspace(0,t,length(U_SingMatMax(1,:)));
D_max = min(length(U_errT(:,1)),31);
plot_chivec = 0:(D_max-1);

%% Plot Entanglement Entropy as a function of time
figure(1)
clf;
subplot(2,1,1)
plot(plot_tvec,U_SingMatMax(1,:),'b','LineWidth',3); hold on;
plot(plot_tvec,V_SingMatMax(1,:),'r','LineWidth',3);
plot(plot_tvec,W_SingMatMax(1,:),'g','LineWidth',3);
axis([0 t 0 max((([U_SingMatMax(1,:),V_SingMatMax(1,:),W_SingMatMax(1,:)]))) ])
title(sprintf('Maximum Entanglement Entropy\nRe = %0.0f, %i x %i x %i grid',Re,2^plot_L,2^plot_L,2^plot_L))
xlabel('t (seconds)');
ylabel('S_{max}(t)');
legend(' V_{x}',' V_{y}', 'V_{z}');


set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')

subplot(2,1,2)
plot(plot_tvec,U_SingMatMax(2,:),'b','LineWidth',3); hold on;
plot(plot_tvec,V_SingMatMax(2,:),'r','LineWidth',3);
plot(plot_tvec,W_SingMatMax(2,:),'g','LineWidth',3);
axis([0 t 1 (plot_L-1) ])
title(sprintf('Bond at which Entanglement Entropy is Maximum of flow field\nRe = %0.1g, %i x %i x %i grid',Re,2^plot_L,2^plot_L,2^plot_L))
xlabel('t (seconds)');
ylabel('bond #');
legend(' V_{x}',' V_{y}', 'V_{z}');

set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')

%% Plot Error in Entanglement Entropy as a function of time
figure(2);
%subplot(3,1,1);
% clf;
% subplot(2,1,1)
% norm = plot_L-1;
% surf(plot_tvec,plot_chivec,U_errT(plot_chivec+1,:)./norm ); view(120,20);
% axis([0 t 0 D_max 0 1 ])
% %title(sprintf('Normalised Cumulative sum of Singular Values squared of V_{x}\n up to each chi, across t. Re = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
% title(sprintf('(a)'));
% xlabel('t (seconds)');
% ylabel('\chi');
% %zlabel('Normalised CumulativeSum(\chi,t)');
% zlabel('Error(\chi,t)')
% 
% set(gca,'FontSize',24,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')
% 
% subplot(2,1,2)
% norm = plot_L-1;
% surf(plot_tvec,plot_chivec,V_errT(plot_chivec+1,:)./norm ); view(120,20);
% axis([0 t 0 D_max 0 1 ])
% %title(sprintf('Normalised Cumulative sum of Singular Values squared of V_{y}\n up to each chi, across t. Re = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
% title(sprintf('(b)'));
% xlabel('t (seconds)');
% ylabel('\chi');
% %zlabel('Normalised CumulativeSum(\chi,t)');
% zlabel('Error(\chi,t)')
% 
% set(gca,'FontSize',24,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')
norm = plot_L-1;
error = ((U_errT(plot_chivec+1,:) + V_errT(plot_chivec+1,:)+W_errT(plot_chivec+1,:) )/(3*norm));
surf(plot_tvec,plot_chivec, error); view(120,20);
ticks = -15:3:0;
%shading interp;

cb = colorbar('Ticks',ticks,'TickLabels',{ num2str(10.^ticks,'%1.1e\n') });
set(cb, 'ylim', [-15 0])

%axis([0 t 0 D_max 0 1 ])
%title(sprintf('Normalised Cumulative sum of Singular Values squared of V_{x}\n up to each chi, across t. Re = %0.1g, %i x %i grid',Re,2^plot_L,2^plot_L))
title(sprintf('(a)'));
xlabel('t (seconds)');
ylabel('\chi');
%zlabel('Normalised CumulativeSum(\chi,t)');
zlabel('log10(Error(\chi,t)')

set(gca,'FontSize',24,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold')