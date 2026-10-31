%calculates the aggregate variables of total kinetic energy, total
%momementum, total vorticity and total divergence. Normalises by dividing these quantities
%by the total number of lattice points.

function [Ek,P,W,div] = calcAggQuantities2D(U,V)

N = (length(U(1,:)));
P = zeros(1,2);
dx = 1/N;

Ek = 0.5*sum(sum(U.^2+V.^2))/N.^2; %Kinetic Energy
P(1) = sum(sum(U))/N.^2; 
P(2) = sum(sum(V))/N.^2;

%Calculate derivatives

dUdx = centrDiff1(U,dx,2);
dVdx = centrDiff1(V,dx,2);
dUdy = centrDiff1(U,dx,1);
dVdy = centrDiff1(V,dx,1);

W = dUdy-dVdx; W = 0.5*sum(sum(W.^2))/N.^2; %Enstrophy
div = dUdx + dVdy; div = sqrt(sum(sum((div.^2)/N^2))); %Divergence RMS

end