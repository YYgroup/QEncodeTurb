%Reconstructs function from MPS in Nik's form

function F = reconstructFromNik3D(MPS)

L = length(MPS);

%Contract the MPS
nodeL = MPS{1};
k=2;
while k ~= length(MPS)+1
    nodeR = MPS{k};
    %If there is a mismatch between the number of columns of nodeL and rows
    %of nodeR, reshape them by bringing the relevant multindex either "up"
    %(into columns) or "down" (into rows). Afterwar the contraction,
    %reverse the operation by bringing the indices of the tensor back to 
    %where they were.
    if (size(nodeL,2)*8 == size(nodeR,1))
        nodeR = rowExpandMat(nodeR',8)';
        nodeL = nodeL*nodeR;
        nodeL = rowExpandMat(nodeL,8);
    elseif (size(nodeL,2) == 8*size(nodeR,1))
        nodeL = rowExpandMat(nodeL,8);
        nodeL = nodeL*nodeR;
        nodeR = rowExpandMat(nodeR,8);
        disp("You should never see this -- something is wrong.");
    else
        nodeL = nodeL*nodeR;
    end
    k = k+1;
end

%Switch the order of the lengthscales to get back
while(size(nodeL,1) ~= 1)
    nodeL = rowExpandMat(nodeL',8)';
end

%C_vec = nodeL;

%MPS is now contracted into "nodeL" and in form: (l_1^x,l_1^y,l_1^z)(l_2^x,l_2^y,l_2^z)...
%(l_L^xl_L^yl_L^z). Now reshape into (l_1^x,l_2^x...l_L^x)(l_1^y,
%l_2^y...l_L^y)(l_1^z,l_2^z...l_L^z)

l_vec = zeros(1,3*L); %Binary number used to loop through all possible
%combinations of l_1^xl_1^yl_1^z...l_L^xl_L^yl_L^z

%Define the indexmap of the current tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(3*L-1));

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = zeros(1,3*L);
F_indexMap(1:3:3*L) = 2.^(L - (1:1:L));
F_indexMap(2:3:3*L) = 2.^(2*L-(1:1:L));
F_indexMap(3:3:3*L) = 2.^(3*L-(1:1:L));

F_vec = zeros(1,8^L);

for k = 1:(8^L)
    F_index = 1 + F_indexMap*l_vec';
    C_index = 1 + C_indexMap*l_vec';
    F_vec(F_index) = nodeL(C_index);
    
    l_vec=binInc(l_vec); %Go through all possible combinations of l by simply doing
    %successive binary increments on it.

end

%Now reshape F_vec into matrix
Nx = 2^L;
Ny = 2^L;
Nz = 2^L;
F= zeros(Nx,Ny,Nz);
for k = 1:Nz
	for j = 1:Ny
        index = 1 + (j-1)*Nx + (k-1)*Nx*Ny;
        F(:,j,k) = F_vec( index:(index+Nx-1) );
    end
end

end