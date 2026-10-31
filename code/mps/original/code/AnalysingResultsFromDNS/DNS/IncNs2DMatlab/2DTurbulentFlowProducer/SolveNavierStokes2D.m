%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% advection-diffusion stage with not(RK4), RK2. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uold = U; Vold = V;
Uc   = U; Vc   = V;
%
Phiold = Phi;
Phic = Phi;

for(rk=1:2)%4)

[dUconv,dVconv] = Compute2DadvectionSkew_MirrorTNT(U,V,D1,KX,KY,dt,dx);
[dUdiff,dVdiff] = Compute2Ddiffusion_MirrorTNT(U,V,D2,KX,KY,Re,dt,dx);
dU = dUconv+dUdiff;
dV = dVconv+dVdiff;
[dU,dV] = project(dU,dV,KX,KY);
%
dPhi = ComputefPhi_PS(U, V, Phi, KX, KY, Pe, dt, dx);


if(rk<2)%4)
    U=Uold+b(rk)*dU;
    V=Vold+b(rk)*dV;
    %
    Phi=Phiold+b(rk)*dPhi;
end

Uc=Uc+a(rk)*dU;
Vc=Vc+a(rk)*dV;
%
Phic=Phic+a(rk)*dPhi;

end

U=Uc; V=Vc;
[U,V] = project(U,V,KX,KY); 
%
Phi=Phic;