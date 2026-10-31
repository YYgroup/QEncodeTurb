function dPhi = ComputefPhi_PS(U, V, Phi, KX, KY, Pe,dt, dx);

% Pe is Peclet number: Re*Pr

dPhidx = centrDiff1(Phi,dx,2);
dPhidy = centrDiff1(Phi,dx,1);

dUPhidx = centrDiff1((U.*Phi),dx,2);
dVPhidy = centrDiff1((V.*Phi),dx,1);


dPhi = 1/Pe*dt*( centrDiff2(Phi,dx,2) + centrDiff2(Phi,dx,1) ) ...
    - dt*.5*( U.*dPhidx + V.*dPhidy + dUPhidx + dVPhidy );


end