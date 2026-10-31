%Reconstructs function from MPS in Nik's form

function F = reconstructFromNik2D(MPS)

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
    if (size(nodeL,2)*4 == size(nodeR,1))
        nodeR = rowExpandMat(nodeR',4)';
        nodeL = nodeL*nodeR;
        nodeL = rowExpandMat(nodeL,4);
    elseif (size(nodeL,2) == 4*size(nodeR,1))
        nodeL = rowExpandMat(nodeL,4);
        nodeL = nodeL*nodeR;
        nodeR = rowExpandMat(nodeR,4);
        disp("You should never see this -- something is wrong.");
    else
        nodeL = nodeL*nodeR;
    end
    k = k+1;
end

%Switch the order of the lengthscales to get back
while(size(nodeL,1) ~= 1)
    nodeL = rowExpandMat(nodeL',4)';
end

%C_vec = nodeL;

%MPS is now contracted into "nodeL" and in form: (l_1^x,l_1^y)(l_2^x,l_2^y)...
%(l_L^xl_L^y). Now reshape into (l_1^x,l_2^x...l_L^x)(l_1^y,
%l_2^y...l_L^y).

l_vec = zeros(1,2*L); %Binary number used to loop through all possible
%combinations of l_1^xl_1^y...l_L^xl_L^y

%Define the indexmap of the current tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(2*L-1));

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = zeros(1,2*L);
F_indexMap(1:2:2*L) = 2.^-(0:1:(L-1));
F_indexMap(2:2:2*L) = 2.^(L-(0:1:(L-1)));
F_indexMap = F_indexMap*2^(L-1);

F_vec = zeros(1,4^L);

for k = 1:(4^L)
    F_index = 1 + F_indexMap*l_vec';
    C_index = 1 + C_indexMap*l_vec';
    F_vec(F_index) = nodeL(C_index); %C_vec = nodeL
    
    l_vec=binInc(l_vec); %Go through all possible combinations of l by simply doing
    %successive binary increments on it.

end

%Now reshape F_vec into matrix
Nx = 2^L;
Ny = 2^L;
F= zeros(Nx,Ny);
for j = 1:Ny
    i = 1 + (j-1)*Nx;
    F(:,j) = F_vec(i:(i-1+Nx));
end

end