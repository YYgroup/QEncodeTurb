%Computes diffusion term

function [dU,dV, dW] = Compute3Ddiffusion_MirrorTNT(U,V,W,D2,KX,KY,KZ,FU, FV, FW, Re,dt,dx)


dU = 1/Re*dt*( centrDiff2(U,dx,2) + centrDiff2(U,dx,1) + centrDiff2(U,dx,3) );
dV = 1/Re*dt*( centrDiff2(V,dx,2) + centrDiff2(V,dx,1) + centrDiff2(V,dx,3) );
dW = 1/Re*dt*( centrDiff2(W,dx,2) + centrDiff2(W,dx,1) + centrDiff2(W,dx,3) );

%dU = 1/Re*dt*(-KX.*KX.*FU + -KY.*KY.*FU + -KZ.*KZ.*FU );
%dV = 1/Re*dt*(-KX.*KX.*FV + -KY.*KY.*FV + -KZ.*KZ.*FV );
%dW = 1/Re*dt*(-KX.*KX.*FW + -KY.*KY.*FW + -KZ.*KZ.*FW );
end