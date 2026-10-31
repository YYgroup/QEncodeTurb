%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for initializing simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L0  = 1; % domain side length
dx = L0/N; % grid spacing assumed dx=dy

%%%%%%%%%%%%%%%%%%%%%%
% create the 3D grid %
%%%%%%%%%%%%%%%%%%%%%%
x = (L0-dx)*(0:(N-1))/(N-1); % to be precise, dx is subtracted here
y = (L0-dx)*(0:(N-1))/(N-1); % otherwise the derivative of e.g. sin(x) wouldnt be differentiable on this periodic grid
z = (L0-dx)*(0:(N-1))/(N-1); %
[X,Y,Z] = meshgrid(x,y,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define starting field     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U  = -U0*sin(k0*X).*cos(k0*Y).*cos(k0*Z);
V =  U0*cos(k0*X).*sin(k0*Y).*cos(k0*Z);
W =  0*X;
%

%V = 0*X;


%% Gaussian profile
% U = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01 + -(Z-0.5).^2/0.01);
% V = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01 + -(Z-0.5).^2/0.01);
% W = exp( -(X-0.5).^2/0.01 + -(Y-0.5).^2/0.01 + -(Z-0.5).^2/0.01);

%[U,V,W] = project3D(U,V,W,KX,KY,KZ,AA);