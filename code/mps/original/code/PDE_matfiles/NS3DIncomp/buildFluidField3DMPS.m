function [MPS_U_TNT, MPS_V_TNT, MPS_W_TNT, U, V, W] = buildFluidField3DMPS(L_orig, chi, k0, U0)
%Builds MPS from either an unstable 3D Fluid Field or Taylor-Green vortex. 
%The Field is built up at a lengthscale corresponding to L_orig, such that
%there are 2^L_orig x 2^L_orig x 2^L_orig grid points.

%Define needed variables
N           = 2^L_orig; % N*N grid, table allows looping over a mesh sensitivity test    

%Create field
CreateFields3D; %uses k0, U0

%Put them into Nikita Gourianov's MPS form and tranpose the U & V
%fields. Transpose because the x directions is defined along the rows in my program, while
%in CreateFields3D x it's defined along cols in accordance with the matlab
%convention.
MPS_U = createMPS3D(permute(U,[2,1,3]),chi);
MPS_V = createMPS3D(permute(V,[2,1,3]),chi);
MPS_W = createMPS3D(permute(W,[2,1,3]),chi);

%Transform into TNT MPS form
MPS_U_TNT = convertMPS_NikTNT(MPS_U,8);
MPS_V_TNT = convertMPS_NikTNT(MPS_V,8);
MPS_W_TNT = convertMPS_NikTNT(MPS_W,8);

end
