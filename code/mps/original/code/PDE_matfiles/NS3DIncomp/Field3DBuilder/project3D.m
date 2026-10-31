function [dU,dV,dW] = project3D(dU,dV,dW,KX,KY,KZ,AA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projects the input field dU,dV,dW (or U,V,V)
% and returns dU,dV pair which is divergence-free.
% The idea is to ensure that the continuity <k, u_k> = 0 
% equation in Fourier transformed form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fft of a field with divergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FU = AA.*fftn(dU);
FV = AA.*fftn(dV);
FW = AA.*fftn(dW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection so that the wave vector  (kx,ky) is orthogonal to (fu,fv,fw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KF = (KX.*FU + KY.*FV + KZ.*FW)./max(1e-10,(KX.^2 + KY.^2 + KZ.^2));

FUU = FU-KF.*KX;
FVV = FV-KF.*KY;
FWW = FW-KF.*KZ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse Fourier transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dU = real(ifftn(FUU));
dV = real(ifftn(FVV));
dW = real(ifftn(FWW));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now (dU,dV,dW) is div free
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
