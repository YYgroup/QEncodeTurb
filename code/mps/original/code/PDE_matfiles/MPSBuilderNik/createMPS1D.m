%Creates a MPS representation of the function F(x), declared as an N
%array. N must be some integer power of 2.

%The maximum bond dimension of the MPS is declared as "chi".

function [MPS,SingMats] = createMPS1D(F,chi)
%% Initialise

%Extract information
N = length(F);
L = log2(N); %Must be integer and equal to log2(N).


%% Reshape into multidimensional array
%F is a multidimensional array indexed as (l_L,l_L-1,...l_1).
%Change indexing to (l_1),(l_2),...(l_L).
F_vec = F;

l_vec = zeros(1,L); %Binary number used to loop through all possible
%combinations of l_1,l_2,...l_L

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(L-1));

%Define the indexmap of the original tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = 2.^(L-(1:L));

%Loop through and define a transformation-matrix T, which maps from the old
%indices to the new one. Then use this to build C. [I chose to do this in
%such an indirect way because then it's easy to map back-and-forth
%between the two index-bases.

C_vec = zeros(1,2^L);

for k = 1:(2^L)
    F_index = 1 + F_indexMap*l_vec';
    C_index = 1 + C_indexMap*l_vec';
    l_vec=binInc(l_vec); %Go through all possible combinations of l by simply doing
    %successive binary increments on it.

    C_vec(C_index) = F_vec(F_index);
    
end

%% Build MPS
%Set the current MPS node to be the initial tensor, represented as C_vec.
%It will be SVD'ed s.t. every node corresponds to one of the "length-
%-scales", ie by one of the multindices.

%Local basis now:
%|l_k> = |l_k>:
%|0> <=> x=+0
%|1> <=> x=+2^-k

node = C_vec;
MPS = {};
SingMats = {};

itNum = 1;
while(itNum < L)
    
    node = rowExpandMat(node,2);
    [U,S,V] = svds(node,chi);
    
    %U is now a proper MPS node -- store into network
    MPS{itNum} = U;
    
    %Save S, then multiply it into V' before continuing.
    SingMats{itNum} = diag(S);
    node = S*V';
    
    itNum = itNum +1;
    
end

node = rowExpandMat(node,2);
MPS{L} = node;

end