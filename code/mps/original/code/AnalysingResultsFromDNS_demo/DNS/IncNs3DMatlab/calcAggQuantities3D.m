%calculates the aggregate variables of total kinetic energy, total
%momementum, total divergence and the curl fields. Normalises by dividing these quantities
%by the total number of lattice points.

function [Ek,P,div, RU, RV, RW] = calcAggQuantities3D(U,V,W)

N = (length(U(1,1,:)));
dx = 1/N;
P = zeros(1,3);

Ek = (U.^2+V.^2+W.^2)/N^3;
P(1) = sum(sum(sum(U)))/N^3;
P(2) = sum(sum(sum(V)))/N^3;
P(3) = sum(sum(sum(W)))/N^3;

%Calculate derivatives
dUdx = centrDiff1(U,dx,2);
dVdx = centrDiff1(V,dx,2);
dWdx = centrDiff1(W,dx,2);
%
dUdy = centrDiff1(U,dx,1);
dVdy = centrDiff1(V,dx,1);
dWdy = centrDiff1(W,dx,1);
%
dUdz = centrDiff1(U,dx,3);
dVdz = centrDiff1(V,dx,3);
dWdz = centrDiff1(W,dx,3);

div = sqrt(sum(sum((  (dUdx+dVdy+dWdz).^2))))/N^3;

RU = dWdy - dVdz;
RV = dUdz - dWdx;
RW = dVdx - dUdy;

end
