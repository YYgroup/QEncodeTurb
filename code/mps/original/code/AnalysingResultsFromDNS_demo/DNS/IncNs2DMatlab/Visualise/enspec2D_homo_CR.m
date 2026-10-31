function loglogPlot = enspec2D_homo_CR(U,V,s,titl,linewidth)
%enspec2D finds the energy spectrum of 2D fluid flow-fields
%across length-scales. The flow fields must be represented on a grid which
%is equispaced along all axes. The minimum length-scale of
%the energy spectrum is defined by the axis with the least number of
%grid points. Assumes the dimensionality is a power of 2.

lengthScale = 1;

FU = fftn(U); FV = fftn(V);
N = length(U(:,1,1));
k = (2*pi/lengthScale)*[0:N/2-1 -N/2+1:0];

[KX, KY] = meshgrid(k,k);

%% Initialise
dk = k(2)-k(1);
Kmax = 2*pi/lengthScale*(N/2-1);

Ek = zeros(1,N/2);
k = zeros(1,N/2); %redefine k

E = FU.*conj(FU) + FV.*conj(FV); %By Parseval's theorem

%% Now integrate energy spectrum
index=1;
Klower = 0; Kupper = 2*dk;

while Kupper <= Kmax
    stepfunc = ( sqrt(KX.^2 + KY.^2) <= Kupper ).*( sqrt(KX.^2 + KY.^2) >= Klower );
    Eatk = E.*stepfunc;
    Eatk = sum(Eatk(:))*2^2.5*dk*dk; % area element = k*dtheta*(2*dk) ~ 2*sqrt(2)*dk*(2*dk)
    
    Ek(index) = Eatk;
    k(index) = Klower+dk;
    
    Klower = Klower + dk;
    Kupper = Kupper + dk;
    index=index+1;
end

%Get rid of the zero-elements
Ek = Ek(1:(index-1));
k = k(1:(index-1));

% Normalise & smooth out
Ek = Ek/sum(Ek);
Ek = smooth(Ek,s);
%Ek(Ek<eps^2)=0;

log10compRatio = zeros(1,length(k));
for i = 1:length(Ek)
    compRatio = (i/length(Ek))^2;
    log10compRatio(i) = log10( compRatio );
end

%logk= log10(k);
%logk(1) = (-1+log10(2*pi));

%% Plot
loglogPlot = plot(log10compRatio,log10(Ek),'LineWidth',linewidth);

%set(gca, 'YScale', 'log')

%xticks([(-1+log10(2*pi)) (0+log10(2*pi)) (1+log10(2*pi)) (2+log10(2*pi)) (3+log10(2*pi)) (4+log10(2*pi))] )
%xticklabels({'0','1','10^1','10^2','10^3','10^4'})


title(titl, 'FontSize',30);
%xlabel('k/2\pi');
%ylabel('E_k');
%ax = gca;
%ax.FontSize = 25;

end