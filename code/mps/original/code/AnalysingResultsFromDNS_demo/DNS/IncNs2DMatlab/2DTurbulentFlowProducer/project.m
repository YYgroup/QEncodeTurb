function [dU,dV] = project(dU,dV,KX,KY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projects the input field dU,dV (or U,V)
% and returns dU,dV pair which is divergence-free.
% The idea is to ensure that the continuity <k, u_k> = 0 
% equation in Fourier transformed form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fft of a field with divergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FU = fft2(dU);
FV = fft2(dV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection so that the wave vector  (kx,ky) is orthogonal to (fu,fv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUU = FU-(KX.*FU + KY.*FV).*KX./max(1e-14,(KX.^2 + KY.^2));
FVV = FV-(KX.*FU + KY.*FV).*KY./max(1e-14,(KX.^2 + KY.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse Fourier transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dU = real(ifft2(FUU));
dV = real(ifft2(FVV));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now (dU,dV) is a div free pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%