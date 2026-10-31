%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for initializing simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L  = 1;%2*pi;    % domain side length
dx = L/N;     % grid spacing assumed dx=dy
Pe = Re;

%%%%%%%%%%%%%%%%%%%%
% create the 2D grid
%%%%%%%%%%%%%%%%%%%%
x = (L-dx)*(0:(N-1))/(N-1); % to be precise, dx is subtracted here 
y = (L-dx)*(0:(N-1))/(N-1); % otherwise the derivative of e.g. sin(x) wouldnt be differentiable on this periodic grid
[X,Y] = meshgrid(x,y);     
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize velocity first to 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(N, N);  % x-velocity
V = zeros(N, N);  % y-velocity 

FU = zeros(N,N);  % tables where Fourier-transforms are stored
FV = zeros(N,N);  % e.g. FU = fft2(U)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wave vectors: see HowToUseFFT2.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kx1 = mod(1/2 + (0:(N-1))/N, 1 ) - 1/2;
ky1 = mod(1/2 + (0:(N-1))/N, 1 ) - 1/2;
kx  = kx1*(2*pi/dx); % wavenumbers
ky  = ky1*(2*pi/dx); % wavenumbers

% these are the wave vector modes one gets from fft2
[KX,KY] = meshgrid(kx,ky); 

% Anti-aliasing filter based on the 2/3*kNyq rule
% AA = double(abs(KX) < (2/3)*max(kx)).*double(abs(KY) < (2/3)*max(ky)); 

% whenever fft2 is called we program AA.*fft2(field) in order to remove the 
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
% Periodic
D1(1,N) = -1; D1(1,2) = +1;
D1(N,N-1) = -1; D1(N,1) = +1;
D2(1,N) = +1; D2(1,2) = +1; D2(1,1) = -2;
D2(N,N-1) = +1; D2(N,1) = +1; D2(N,N) = -2;

%Fill elements corresponding to the boundary
% Periodic-reflective
% D1(1,N) = +1; D1(1,2) = +1;
% D1(N,N-1) = -1; D1(N,1) = -1;
% D2(1,N) = -1; D2(1,2) = +1; D2(1,1) = -2;
% D2(N,N-1) = +1; D2(N,1) = -1; D2(N,N) = -2;

%Scale appropriately
D1 = sparse(D1)/(2*dx);
D2 = sparse(D2)/(dx^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta 2/4 coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RK 4
%a1= 1/6; a2 = 1/3; a3 = 1/3; a4 = 1/6;
%b1= 0.5; b2=0.5; b3=1; b4=1;
% RK 2
a1 = 0; a2 = 1; a3 = nan; a4 = nan;
b1 = 0.5; b2 = nan; b3 = nan; b4 = nan;

%disp('warning:euler!');
%a1=1;b1=0; a2=0;

a=[a1 a2 a3 a4]; b=[b1 b2 b3 b4];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Invent a non-linear function to drive the flow into turbulence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % % Non-linear vortices
% U  = -Umax*cos(kmax*X).*sin(kmax*Y);  
% V =  Umax*sin(kmax*X).*cos(kmax*Y);
% 
% U  = -Umax*cos(kmax*X.*V).*sin(kmax*Y.*U);  
% V = Umax*sin(kmax*X.*U).*cos(kmax*Y.*V);
% 
% %[U,V] = project(U,V,KX,KY,AA);

%% Jet
%A1 = 1/max(2*kmax,19/8*alpha);
A1= 0.025;
D = 2/10;

% up and vp are velocity perturbations
up =    2.0.*kmax^2.*(Y -0.5-D/2).*exp( -kmax^2.*(Y-0.5-D/2).^2 ).*( sin(8*pi*alpha*X)+sin(24*pi*alpha*X)+sin(6*pi*alpha*X) );
% up =up+ 2.0.*kmax^2.*(Y -0.5+D/2).*exp( -kmax^2.*(Y-0.5+D/2).^2 ).*( sin(32*pi*alpha*X+ 3*pi/2)+sin(4*pi*alpha*X+ pi/2)+sin(2*pi*alpha*X) );
up =up+ 2.0.*kmax^2.*(Y -0.5+D/2).*exp( -kmax^2.*(Y-0.5+D/2).^2 ).*( sin(8*pi*alpha*X)+sin(24*pi*alpha*X)+sin(6*pi*alpha*X) );

vp =    exp( -kmax^2.*(Y-0.5-D/2).^2 ).*alpha.*( 8*pi*cos(8*pi*alpha.*X)+24*pi.*cos(24*pi*alpha.*X)+6*pi*cos(6*pi*alpha.*X) );
vp =vp+ exp( -kmax^2.*(Y-0.5+D/2).^2 ).*alpha.*( 8*pi*cos(8*pi*alpha.*X)+24*pi.*cos(24*pi*alpha.*X)+6*pi*cos(6*pi*alpha.*X) );
% vp =vp+ exp( -kmax^2.*(Y-0.5+D/2).^2 ).*alpha.*( 32*pi*cos(32*pi.*alpha.*X+ 3*pi/2)+4*pi.*cos(4*pi.*alpha.*X+ pi/2)+2*pi.*cos(2*pi.*alpha.*X) );

U = 1*( tanh(kmax.*(Y-0.5+D/2)) - tanh(kmax.*(Y-0.5-D/2)) ) -1;

pscale = max(max( sqrt(up.^2+vp.^2)));

up = up*(A1/pscale);
vp = vp*(A1/pscale);

U = 0.5*Umax.*(U + up);
V = 0.5*Umax.*(0.*X  + vp);

%Phi:
kmax_phi = 20;
Phi = 1*( tanh(kmax_phi.*(Y-0.5+D/2)) - tanh(kmax_phi.*(Y-0.5-D/2)) );
Phi = Phi./max(max(Phi));

%% Gaussian profile
%U = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01);
%V = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01);

%[U,V] = project(U,V,KX,KY,AA);