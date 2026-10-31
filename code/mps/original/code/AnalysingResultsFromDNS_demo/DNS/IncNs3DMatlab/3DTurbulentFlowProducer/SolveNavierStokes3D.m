%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% advection-diffusion stage with RK4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uold = U; Vold = V; Wold = W;
Uc   = U; Vc   = V; Wc = W; 

for(rk=1:2)%4)

FU = fftn(U); FV = fftn(V); FW = fftn(W);
[dUconv,dVconv, dWconv] = Compute3DadvectionSkew_MirrorTNT(U,V,W,D1,KX,KY,KZ,FU, FV, FW, AA,dt,dx);
[dUdiff,dVdiff, dWdiff] = Compute3Ddiffusion_MirrorTNT(U,V,W,D2,KX,KY,KZ,FU, FV, FW,Re,dt,dx);

dU = dUconv+dUdiff;
dV = dVconv+dVdiff;
dW = dWconv+dWdiff;

[dU,dV,dW] = project3D(dU,dV,dW,KX,KY,KZ,AA); 

if(rk<2)
    U=Uold+b(rk)*dU;
    V=Vold+b(rk)*dV; 
    W=Wold+b(rk)*dW; 
end
    Uc=Uc+a(rk)*dU;
    Vc=Vc+a(rk)*dV;
    Wc=Wc+a(rk)*dW;
end

U=Uc; V=Vc; W=Wc;
[U,V,W] = project3D(U,V,W,KX,KY,KZ,AA); 
