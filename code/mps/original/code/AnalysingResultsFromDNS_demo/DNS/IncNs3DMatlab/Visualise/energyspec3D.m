function [K,EK] = energyspec3D(FU,FV,FW,KX,KY,KZ)
%ENSPEC3D finds the energy spectrum of 3D fluid flow-fields
%across length-scales. The flow fields must be represented on a grid which
%is equispaced along all axes. The minimum length-scale of
%the energy spectrum is defined by the axis with the least number of
%grid points. Assumes the dimensionality is a power of 2.

%% Initialise
dk = KX(1,2,1) - KX(1,1,1);
Kmax = max(abs(KX(:)));

EK = zeros(1,Kmax/dk -1);
K =  zeros(1,Kmax/dk -1);

E = FU.*conj(FU) + FV.*conj(FV) + FW.*conj(FW); %By Parseval's theorem

%% Now integrate energy spectrum
index=1;
Klower = 0; Kupper = dk;

while Kupper <= Kmax
    stepfunc = ( sqrt(KX.^2 + KY.^2 + KZ.^2) <= Kupper ).*( sqrt(KX.^2 + KY.^2 + KZ.^2) >= Klower );
    Eatk = E.*stepfunc;
    
    EK(index) = sum(Eatk(:));
    K(index) = Klower+0.5*dk;
    
    Klower = Klower + dk;
    Kupper = Kupper + dk;
    index=index+1;
end

% Normalise
EK = EK/sum(EK(:));
EK(EK<=eps) = 0;

end