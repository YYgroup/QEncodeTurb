%Calculates the aggregate variables along the column dimensions of the field (U,V). 
%Normalises by dividing these quantities by the number of columns.

function [P, K, O, W, div, E, VortThick, dP1dy,ReStress11,ReStress12, ReStress22] = calcAggAlongColsQuantities2D(U,V)

N = (length(U(1,:)));
dx = 1/N;

% x-avged Momenta
P(:,1) = sum(U,2)/N; 
P(:,2) = sum(V,2)/N;

% x-avged kin energy and turb kin energy 
E = 0.5*sum(U.^2 + V.^2, 2)/N;
K = 0.5*sum((U-P(:,1)).^2 + (V-P(:,2)).^2, 2)/N;

%Calculate derivatives
dUdx = centrDiff1(U,dx,2);
dVdx = centrDiff1(V,dx,2);
dUdy = centrDiff1(U,dx,1);
dVdy = centrDiff1(V,dx,1);
dVtdx = centrDiff1(V-P(:,2),dx,2);
dUtdy = centrDiff1(U-P(:,1),dx,1);

fieldCurl = dUdy-dVdx; 
fieldCurlt = dUtdy-dVtdx; 

% x-avged vorticity functions
O = sum(fieldCurl, 2)/N;

W = sum((fieldCurlt).^2,2)/N;

% x,y-avged divergence
div = dUdx + dVdy; div = sqrt(sum(sum((div.^2)/N^2))); %Divergence RMS

% vorticity thickness and x-avged x momenta y-derivative
DU = -(P(1,1) - P(round(N/2)+1,1));
dP1dy = centrDiff1(P(:,1),dx,1);
VortThick = DU/max(dP1dy);

% x-avged Reynolds stress

ReStress11 = sum( (U-P(:,1)).*(U-P(:,1)) ,2)/N;
ReStress12 = sum( (U-P(:,1)).*(V-P(:,2)) ,2)/N;
ReStress22 = sum( (V-P(:,2)).*(V-P(:,2)) ,2)/N;

end