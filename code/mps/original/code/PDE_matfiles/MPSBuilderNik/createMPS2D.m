%Creates a MPS representation of the function F(x,y), declared as an Nx x Ny
%matrix. (x,y) must be the elements of a rectangular grid with an even 
%number of points along each axis. The x-axis points "downwards" and is 
%proportional to the rows, while the y-axis points towards the right. 
%Nx and Ny must be equal to eachother and be some integer power of 2.

%The maximum bond dimension of the MPS is declared as "chi".

function MPS = createMPS2D(F,chi)
%% Initialise

%Extract information
[Nx,Ny] = size(F); %Assume rows correspond to x, columns to y. Nx must be
%equal to Ny.
L = log2(Nx); %Must be integer and equal to log2(Ny).

%Vectorise function. This function will now be a multidimensional array
%indexed as (l_L^x,l_L-1^x...l_1^x)(l_L^y,l_L-1^y...l_1^y)
F_vec = zeros(1,4^L);
for j = 1:Ny
    i = 1 + (j-1)*Nx;
    F_vec( i:(i-1+Nx) ) = F(:,j);
end

%% Reshape into multidimensional array
%Change indexing from (l_L^x,l_L-1^x...l_1^x)(l_L^y,l_L-1^y...l_1^y) to
%(l_1^x,l_1^y)(l_2^x,l_2^y)...(l_L^xl_L^y).

l_vec = zeros(1,2*L); %Binary number used to loop through all possible
%combinations of l_1^xl_1^y...l_L^xl_L^y

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(2*L-1));

%Define the indexmap of the original tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = zeros(1,2*L);
F_indexMap(1:2:2*L) = 2.^(L-(1:1:L));
F_indexMap(2:2:2*L) = 2.^(2*L-(1:1:L));

%Loop through and define a transformation-matrix T, which maps from the old
%indices to the new one. Then use this to build C. [I chose to do this in
%such an indirect way because then it's easy to map back-and-forth
%between the two index-bases.

C_vec = zeros(1,4^L);

for k = 1:(4^L)
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

%Switch from indexing x and y seperately to combining them into one single
%multindex. Because of this, the hilbert space of local basis is now 4
%instead of 2. Formal definition of new multindex:
%|l_k> = |l_k^x,l_k^y>:
%|0> = |00> <=> x=+0,y=+0
%|1> = |10> <=> x=+2^-k,y=+0
%|2> = |01> <=> x=+0,y=+2^-k
%|3> = |11> <=> x=+2^-k,y=+2^-k

node = C_vec;
MPS = {};
SingMats = {};

itNum = 1;
while(itNum < L)
    
    node = rowExpandMat(node,4);
    [U,S,V] = svds(node,chi);
    
    %U is now a proper MPS node -- store into network
    MPS{itNum} = U;
    
    %Save S, then multiply it into V' before continuing.
    SingMats{itNum} = S;
    node = S*V';
    
    itNum = itNum +1;
    
end

node = rowExpandMat(node,4);
MPS{L} = node;

end