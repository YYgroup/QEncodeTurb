%calculates the aggregate variables of total kinetic energy, total
%momementum, total divergence and the curl fields. Normalises by dividing these quantities
%by the total number of lattice points.

function [Ek,P,div, lambda2, RU, RV, RW] = calcAggQuantities3D(U,V,W)

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

%Divergence
div = sqrt(sum(sum((  (dUdx+dVdy+dWdz).^2))))/N^3;

%Vorticity
RU = dVdx - dUdy;
RV = dUdz - dWdx;
RW = dVdx - dUdy;

% lambda2-criterion
lambda2 = zeros(N^3,1);

%Reshape
dUdx = reshape(dUdx,[N^3,1]);
dUdy = reshape(dUdy,[N^3,1]);
dUdz = reshape(dUdz,[N^3,1]);
dVdx = reshape(dVdx,[N^3,1]);
dVdy = reshape(dVdy,[N^3,1]);
dVdz = reshape(dVdz,[N^3,1]);
dWdx = reshape(dWdx,[N^3,1]);
dWdy = reshape(dWdy,[N^3,1]);
dWdz = reshape(dWdz,[N^3,1]);

for i = 1:N^3
    dudx = dUdx(i); dvdx = dVdx(i); dwdx = dWdx(i);
    dudy = dUdy(i); dvdy = dVdy(i); dwdy = dWdy(i);
    dudz = dUdz(i); dvdz = dVdz(i); dwdz = dWdz(i);
    %
    J = [dudx, dudy, dudz; dvdx, dvdy, dvdz; dwdx, dwdy, dwdz;];
    %
    S = 0.5*(J+J'); 
    Omega = 0.5*(J-J');
    lambdas = sort(eig(S*S + Omega*Omega),'descend');
    %
    lambda2(i) = -lambdas(2);
end

%Reshape lamda2 back
lambda2 = reshape(lambda2,[N,N,N]);

end
