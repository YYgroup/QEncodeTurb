%Reconstructs function from MPS in Nik's form

function F_vec = reconstructFromNik1D(MPS)

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
    if (size(nodeL,2)*2 == size(nodeR,1))
        nodeR = rowExpandMat(nodeR',2)';
        nodeL = nodeL*nodeR;
        nodeL = rowExpandMat(nodeL,2);
    elseif (size(nodeL,2) == 2*size(nodeR,1))
        nodeL = rowExpandMat(nodeL,2);
        nodeL = nodeL*nodeR;
        nodeR = rowExpandMat(nodeR,2);
        disp("You should never see this -- something is wrong.");
    else
        nodeL = nodeL*nodeR;
    end
    k = k+1;
end

%Switch the order of the lengthscales to get back
while(size(nodeL,1) ~= 1)
    nodeL = rowExpandMat(nodeL',2)';
end

%MPS is now contracted into "nodeL" and is in form: (l_1)(l_2)...
%(l_L). Now reshape into (l_L,l_L-1...l_1).

l_vec = zeros(1,L); %Binary number used to loop through all possible
%combinations of l_1,l_2,l_L

%Define the indexmap of the current tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(L-1));

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = 2.^(L-(1:1:L));

F_vec = zeros(2^L,1);

for k = 1:(2^L)
    F_index = 1 + F_indexMap*l_vec';
    C_index = 1 + C_indexMap*l_vec';
    F_vec(F_index) = nodeL(C_index); %C_vec = nodeL
    
    l_vec=binInc(l_vec); %Go through all possible combinations of l by simply doing
    %successive binary increments on it.

end

%Done.

end