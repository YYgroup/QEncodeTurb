%Calculates the aggregate variables along the column dimensions of the field (U,V). 
%Normalises by dividing these quantities by the number of columns.

function [reStress, scalarThick, scalarFlux] = calcAggAlongColsConsScQuantities2D(Phi, U, V)

N = (length(U(1,:)));
dx = 1/N;

% x-avged Momenta & phi
P(:,1) = sum(U,2)/N; 
P(:,2) = sum(V,2)/N;
Phi_avg = sum(Phi,2)/N;

% x-avged Reynolds stress
reStress = sum( (U-P(:,1)).*(V-P(:,2)) ,2)/N;

% scalar thickness and x-avged x momenta y-derivative
DPhi_avg = -( Phi_avg(1) - Phi_avg(round(N/2)+1) );
dPhi_avgdy = centrDiff1(Phi_avg,dx,1);
scalarThick = DPhi_avg/max(dPhi_avgdy);

% x-avged Scalar flux
scalarFlux = sum( (Phi-Phi_avg).*(V-P(:,2)) ,2)/N;

end