%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for initializing simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L  = 2*pi; % domain side length
dx = L/N; % grid spacing assumed dx=dy

%%%%%%%%%%%%%%%%%%%%%%
% create the 3D grid %
%%%%%%%%%%%%%%%%%%%%%%
x = (L-dx)*(0:(N-1))/(N-1); % to be precise, dx is subtracted here
y = (L-dx)*(0:(N-1))/(N-1); % otherwise the derivative of e.g. sin(x) wouldnt be differentiable on this periodic grid
z = (L-dx)*(0:(N-1))/(N-1); %
[X,Y,Z] = meshgrid(x,y,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize velocity first to 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FU = zeros(N,N);  % tables where Fourier-transforms are stored
FV = zeros(N,N);  % e.g. FU = fft2(U)
FW = zeros(N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wave vectors: see HowToUseFFT2.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kx1 = mod(1/2 + (0:(N-1))/N, 1 ) - 1/2;
ky1 = kx1;
kz1 = kx1;

kx  = kx1*(2*pi/dx); % wavenumbers
ky  = ky1*(2*pi/dx); % wavenumbers
kz  = kz1*(2*pi/dx); % wavenumbers

% these are the wave vector modes one gets from fft2
[KX,KY,KZ] = meshgrid(kx,ky,kz);

% Anti-aliasing filter based on the 2/3*kNyq rule
AA = double(abs(KX) < (2/3)*max(kx)).*double(abs(KY) < (2/3)*max(ky)).*double(abs(KZ) < (2/3)*max(kz));

% whenever fft3 is called we multiply with AA, i.e. run as AA.*fft3(field), in order to remove the
% aliasing error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First and second order differencing operators %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise
D1 = zeros(N,N);
D2 = D1;

%Fill the bulk of the matrix
for i = 2:(N-1)
    D1(i,i-1) = -1;
    D1(i,i+1) = +1;
    
    D2(i,i)   = -2;
    D2(i,i-1) = +1;
    D2(i,i+1) = +1;
end


%Fill elements corresponding to the boundary
D1(1,N) = -1; D1(1,2) = +1;
D1(N,N-1) = -1; D1(N,1) = +1;
D2(1,N) = +1; D2(1,2) = +1; D2(1,1) = -2;
D2(N,N-1) = +1; D2(N,1) = +1; D2(N,N) = -2;

%Scale appropriately
D1 = D1/(2*dx);
D2 = D2/(dx^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta !4, 2 actually, coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%a1= 1/6; a2 = 1/3; a3 = 1/3; a4 = 1/6;
%b1= 0.5; b2=0.5; b3=1; b4=1;
a1 = 0; a2 = 1; a3 = nan; a4 = nan;
b1 = 0.5; b2 = nan; b3 = nan; b4 = nan;
a=[a1 a2 a3 a4]; b=[b1 b2 b3 b4];

%% TG Vortex
U  = -Umax*sin(kmax*X).*cos(kmax*Y).*cos(kmax*Z);
V =  Umax*cos(kmax*X).*sin(kmax*Y).*cos(kmax*Z);
W =  0*X;
%

%% Invent a non-linear function to drive the flow into turbulence
% 
% U  = -Umax*cos(kmax*X).*sin(kmax*Y).*sin(kmax*Z);
% V =  0.5*Umax*sin(kmax*X).*cos(kmax*Y).*sin(kmax*Z);
% W =  0.5*Umax*sin(kmax*X).*sin(kmax*Y).*cos(kmax*Z);
% 
% U  = -Umax*cos(kmax*X.*V).*sin(kmax*Y.*W).*sin(kmax*Z.*U);
% V = Umax*sin(kmax*X.*W).*cos(kmax*Y.*U).*sin(kmax*Z.*V);
% W =  Umax*sin(kmax*X.*U).*cos(kmax*Y.*V).*cos(kmax*Z.*W);
% 
% [U,V,W] = project3D(U,V,W,KX,KY,KZ,AA);

%% Gaussian profile
%  U = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01 + -(Z-0.5).^2/0.01);
%  V = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01 + -(Z-0.5).^2/0.01);
%  W = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01 + -(Z-0.5).^2/0.01);

