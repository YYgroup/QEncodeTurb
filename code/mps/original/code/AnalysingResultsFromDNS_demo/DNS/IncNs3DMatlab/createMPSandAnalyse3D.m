%Creates a MPS representation of the function F(x,y,z), declared as an Nx x
%Ny x Nz matrix. (x,y,z) must be the elements of a rectangular grid with an even 
%number of points along each axis. The x-axis points "downwards" and is 
%proportional to the rows, while the y-axis points towards the right. 
%Nx and Ny must be equal to eachother and be some integer power of 2.

%The maximum bond dimension of the MPS is declared as "chi".

function [SingVals,entanglementEntropies] = createMPSandAnalyse3D(F,chi)
%% Initialise

%Extract information
[Nx,Ny,Nz] = size(F); %Assume rows correspond to x, columns to y, and vert
%-ical axis to z. Assume Nx = Ny = Nz.
L = log2(Nx); %Must be integer and equal to log2(Ny).

%Vectorise function. This function will now be a multidimensional array
%indexed as (l_L^x,l_L-1^x...l_1^x)(l_L^y,l_L-1^y...l_1^y)(l_L^z,l_L-1^z...l_1^z)
F_vec = zeros(1,8^L);
for k = 1:Nz
    for j = 1:Ny
        index = 1 + (j-1)*Nx + (k-1)*Nx*Ny;
        F_vec( index:(index+Nx-1) ) = F(:,j,k);
    end
end

%% Reshape into multidimensional array
%Change indexing from (l_L^x,l_L-1^x...l_1^x)(l_L^y,l_L-1^y...l_1^y)(l_L^z,l_L-1^z...l_1^z) to
%(l_1^x,l_1^y,l_1^z)(l_2^x,l_2^y,l_2^z)...(l_L^x,l_L^y,l_L^z).

l = zeros(1,3*L); %Binary number used to loop through all possible
%combinations of l_1^xl_1^yl_1^<...l_L^xl_L^yl_L^z

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(3*L-1));

%Define the indexmap of the original tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = zeros(1,3*L);
F_indexMap(1:3:3*L) = 2.^(L - (1:1:L));
F_indexMap(2:3:3*L) = 2.^(2*L-(1:1:L));
F_indexMap(3:3:3*L) = 2.^(3*L-(1:1:L));

%Loop through and define a transformation-matrix T, which maps from the old
%indices to the new one. Then use this to build C. [I chose to do this in
%such an indirect way because then it's easy to map back-and-forth
%between the two index-bases.

C_vec = zeros(1,8^L);

for k = 1:(8^L)
    F_index = 1 + F_indexMap*l';
    C_index = 1 + C_indexMap*l';
    
    C_vec(C_index) = F_vec(F_index);
    
    l=binInc(l); %Go through all possible combinations of l by simply doing
    %successive binary increments on it.
end

%% Build MPS
%Set the current MPS node to be the initial tensor, represented as C_vec.
%It will be SVD'ed s.t. every node corresponds to one of the "length-
%-scales", ie by one of the multindices.
node = C_vec;
MPS = {};
SingVals = zeros(L-1,chi);
entanglementEntropies = zeros(2,L-1);

itNum = 1;
while(itNum ~= L)
    
    node = rowExpandMat(node,8);
    [U,S,V] = svds(node,chi);
        
    %U is now a proper MPS node -- store into network
    MPS{itNum} = U;
    
    %Save S, then multiply it into V' before continuing.
    node = S*V';
        
    %Save singular values
    SingVals(itNum,1:length(diag(S))) = diag(S);
    
    %Calculate the entanglement entropy at each level, along with its theoretical maximum value (for scale).
    p = S^2./trace(S^2);
    entEnt = -diag(p)'*log2(diag(p));
    entanglementEntropies(1,itNum) = entEnt;
    entanglementEntropies(2,itNum) = log2(size(p,1));
    
    itNum = itNum +1;
    
end

node = rowExpandMat(node,8);
MPS{L} = node;
