%Computes advection term

function [dU,dV,dW] = Compute3DadvectionSkew_MirrorTNT(U,V,W,D1,KX,KY,KZ,FU, FV, FW, AA,dt, dx)

dUdx = centrDiff1(U,dx,2); dUdy = centrDiff1(U,dx,1); dUdz = centrDiff1(U,dx,3);
dVdx = centrDiff1(V,dx,2); dVdy = centrDiff1(V,dx,1); dVdz = centrDiff1(V,dx,3);
dWdx = centrDiff1(W,dx,2); dWdy = centrDiff1(W,dx,1); dWdz = centrDiff1(W,dx,3);

 
dUUdx = centrDiff1((U.*U),dx,2);
dUVdx = centrDiff1(U.*V,dx,2);
dUWdx = centrDiff1(U.*W,dx,2);
dVUdy = centrDiff1(V.*U,dx,1);
dVVdy = centrDiff1(V.*V,dx,1);
dVWdy = centrDiff1(V.*W,dx,1);
dWUdz = centrDiff1(W.*U,dx,3);
dWVdz = centrDiff1(W.*V,dx,3);
dWWdz = centrDiff1(W.*W,dx,3);


dU = -dt/2*(U.*dUdx + V.*dUdy + W.*dUdz + dUUdx + dVUdy + dWUdz);
dV = -dt/2*(U.*dVdx + V.*dVdy + W.*dVdz + dUVdx + dVVdy + dWVdz);
dW = -dt/2*(U.*dWdx + V.*dWdy + W.*dWdz + dUWdx + dVWdy + dWWdz);


%dU = (-U.*real(ifftn(i*KX.*AA.*FU)) - V.*real(ifftn(i*KY.*AA.*FU)) -W.*real(ifftn(i*KZ.*AA.*FU)) - ...
%   real(ifftn(i*KX.*AA.*fftn(U.*U))) - real(ifftn(i*KY.*AA.*fftn(U.*V))) - real(ifftn(i*KZ.*AA.*fftn(U.*W))) )*dt/2;
%dV = (-U.*real(ifftn(i*KX.*AA.*FV)) - V.*real(ifftn(i*KY.*AA.*FV)) -W.*real(ifftn(i*KZ.*AA.*FV)) - ...
%   real(ifftn(i*KX.*AA.*fftn(V.*U))) - real(ifftn(i*KY.*AA.*fftn(V.*V))) - real(ifftn(i*KZ.*AA.*fftn(V.*W))) )*dt/2;
%dW = (-U.*real(ifftn(i*KX.*AA.*FW)) - V.*real(ifftn(i*KY.*AA.*FW)) -W.*real(ifftn(i*KZ.*AA.*FW)) - ...
%   real(ifftn(i*KX.*AA.*fftn(W.*U))) - real(ifftn(i*KY.*AA.*fftn(W.*V))) - real(ifftn(i*KZ.*AA.*fftn(W.*W))) )*dt/2;

end