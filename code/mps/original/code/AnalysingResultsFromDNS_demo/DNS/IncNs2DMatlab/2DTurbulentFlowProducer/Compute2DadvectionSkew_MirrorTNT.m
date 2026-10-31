%Computes advection term

function [dU,dV] = Compute2DadvectionSkew_MirrorTNT(U,V,D1,KX,KY,dt,dx)

dUdx = centrDiff1(U,dx,2); dUdy = centrDiff1(U,dx,1);
dVdx = centrDiff1(V,dx,2); dVdy = centrDiff1(V,dx,1);

 
dUUdx = centrDiff1((U.*U),dx,2);
dUUdy = centrDiff1((U.*U),dx,1);
dUVdx = centrDiff1(U.*V,dx,2);
dUVdy = centrDiff1(U.*V,dx,1);
dVVdx = centrDiff1(V.*V,dx,2); 
dVVdy = centrDiff1(V.*V,dx,1);


%dU = -dt*( dUUdx + dUVdy);
%dV = -dt*( dUVdx + dVVdy);

 dU = -dt/2*(U.*dUdx + V.*dUdy + dUUdx + dUVdy);
 dV = -dt/2*(U.*dVdx + V.*dVdy + dUVdx + dVVdy);

%dU = ( -U.*real(ifft2(i*KX.*fft2(U))) - V.*real(ifft2(i*KY.*fft2(U))) - ...
%      real(ifft2(i*KX.*fft2(U.*U))) - real(ifft2(i*KY.*fft2(U.*V))) )*dt/2;
%dV = ( -U.*real(ifft2(i*KX.*fft2(V))) - V.*real(ifft2(i*KY.*fft2(V))) - ...
%      real(ifft2(i*KX.*fft2(U.*V))) - real(ifft2(i*KY.*fft2(V.*V))) )*dt/2;

end