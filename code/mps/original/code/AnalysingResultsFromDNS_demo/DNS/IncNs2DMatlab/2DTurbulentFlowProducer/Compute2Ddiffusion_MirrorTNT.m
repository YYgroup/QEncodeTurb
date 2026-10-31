%Computes diffusion term

function [dU,dV] = Compute2Ddiffusion_MirrorTNT(U,V,D2,KX,KY,Re,dt,dx)

%dU = 1/Re*dt*( real(ifft2( -KX.*KX.*fft2(U))) + real(ifft2( -KY.*KY.*fft2(U))));
%dV = 1/Re*dt*( real(ifft2( -KX.*KX.*fft2(V))) + real(ifft2( -KY.*KY.*fft2(V))));

dU = 1/Re*dt*( centrDiff2(U,dx,2) + centrDiff2(U,dx,1) );
dV = 1/Re*dt*( centrDiff2(V,dx,2) + centrDiff2(V,dx,1) );


end