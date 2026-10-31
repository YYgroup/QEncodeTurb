function [MPS_U_TNT, MPS_V_TNT, U, V] = buildFluidField2DMPS(Umax, L_orig, chi, kmax, alpha)
%Builds MPS from an unstable 2D Fluid Field. The Field is built up at a
%lengthscale corresponding to L_orig (# grid points = 2^L_orig x 2^L_orig).

%Define needed variables
N = 2^L_orig; % N*N grid, table allows looping over a mesh sensitivity test    

%Make length is non-dimensionalised
L = 1;

%Create field
CreateFields2D; %uses L

% %Rescale into non-dimensional units
% U = U/Umax;
% V = V/Umax;

%Put them into Nikita Gourianov's MPS form (and tranpose the U & V
%fields as the x directions is defined along the rows in my program, while
%in CreateFields2D they define x along cols)
MPS_U = createMPS2D(U',chi);
MPS_V = createMPS2D(V',chi);

%Transform into TNT MPS form
MPS_U_TNT = convertMPS_NikTNT(MPS_U,4);
MPS_V_TNT = convertMPS_NikTNT(MPS_V,4);

end