%Creates a MPS representation of the function F(x,y), declared as an Nx x Ny
%matrix. (x,y) must be the elements of a rectangular grid with an even 
%number of points along each axis. The x-axis points "downwards" and is 
%proportional to the rows, while the y-axis points towards the right. 
%Nx and Ny must be equal to eachother and be some integer power of 2.

%The maximum bond dimension of the MPS is declared as "chi".

function [SingMats, entanglementEntropies] = createMPS(F,chi)
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

l = zeros(1,2*L); %Binary number used to loop through all possible
%combinations of l_1^xl_1^y...l_L^xl_L^y

%Define the indexmap of the reshaped tensor. When it's multiplied with l_vec,
%it gives the index of C corresponding with the combination of l_vec (where
%indexing starts at 0).
C_indexMap = 2.^(0:1:(2*L-1));

%Define the indexmap of the original tensor. When it's multiplied with l_vec,
%it gives the index of F corresponding with the combination of l_vec (where
%indexing starts at 0).
F_indexMap = zeros(1,2*L);
F_indexMap(1:2:2*L) = 2.^-(0:1:(L-1));
F_indexMap(2:2:2*L) = 2.^(L-(0:1:(L-1)));
F_indexMap = F_indexMap*2^(L-1);

%Loop through and define a transformation-matrix T, which maps from the old
%indices to the new one. Then use this to build C. [I chose to do this in
%such an indirect way because then it's easy to map back-and-forth
%between the two index-bases.

%Option 1 only
C_vec = zeros(1,4^L);

%Option 2 only
%T = spalloc(4^L,4^L,4^L);

for k = 1:(4^L)
    F_index = 1 + F_indexMap*l';
    C_index = 1 + C_indexMap*l';
    l=binInc(l); %Go through all possible combinations of l by simply doing
    %successive binary increments on it.

    %Option 1:
    C_vec(C_index) = F_vec(F_index);

    %Option 2
    % T(C_index,F_index) = 1;
    % if(mod(k,10^4)==0)
    %     fprintf("%.1f Percent Ready\n",100*k/(4^L));
    % end
end

%Option 2 only
%C_vec = (T*F_vec')';

%Switch from indexing x and y seperately to combining them into one single
%multindex. Because of this, the hilbert space of local basis is now 4
%instead of 2. Formal definition of new multindex:
%|l_k> = |l_k^x,l_k^y>:
%|0> = |00> <=> Dx=+0,Dy=+0
%|1> = |10> <=> Dx=+2^-k,Dy=+0
%|2> = |01> <=> Dx=+0,Dy=+2^-k
%|3> = |11> <=> Dx=+2^-k,Dy=+2^-k


%% Build MPS
%Set the current MPS node to be the initial tensor, represented as C_vec.
%It will be SVD'ed s.t. every node corresponds to one of the "length-
%-scales", ie by one of the multindices.
node = C_vec;
MPS = {};
SingMats = {};

itNum = 1;
while(itNum ~= L)
    
    node = rowExpandMat(node,4);
    [U,S,V] = svds(node,chi);
    
    %U is now a proper MPS node -- store into network
    MPS{itNum} = U;
    
    %Save S, then multiply it into V' before continuing.
    SingMats{itNum} = S;
    node = S*V';
    
    itNum = itNum +1;
    
end

MPS{L} = node;

%% Calculate Entanglement Entropy
%Calculate the entanglement entropy at each level, along with its
%theoretical maximum value (for scale).
entanglementEntropies = zeros(2,L-1);

%First normalise the content of SingMats Squared, to get a QM-style
%probabilistic interpretation. Then calculate the entanglement entropy.
for itNum = 1:L-1
    S = SingMats{itNum};
    %All calculations can be done through simple matrix multiplications
    %due to S being diagonal.
    S = S^2./trace(S^2);
    entEnt = -diag(S)'*log2(diag(S));
    entanglementEntropies(1,itNum) = entEnt;
    entanglementEntropies(2,itNum) = log2(size(S,1));
end

%% Reconstruct function from MPS to ensure everything works correctly
% 
% nodeL = MPS{1};
% k=2;
% while k ~= length(MPS)+1
%     nodeR = MPS{k};
%     %If there is a mismatch between the number of columns of nodeL and rows
%     %of nodeR, reshape them by bringing the relevant multindex either "up"
%     %(into columns) or "down" (into rows). Afterwar the contraction,
%     %reverse the operation by bringing the indices of the tensor back to
%     %where they were.
%     if (size(nodeL,2)*4 == size(nodeR,1))
%         nodeR = rowExpandMat(nodeR',4)';
%         nodeL = nodeL*nodeR;
%         nodeL = rowExpandMat(nodeL,4);
%     elseif (size(nodeL,2) == 4*size(nodeR,1))
%         nodeL = rowExpandMat(nodeL,4);
%         nodeL = nodeL*nodeR;
%         nodeR = rowExpandMat(nodeR,4);
%         disp("hei");
%     else
%         nodeL = nodeL*nodeR;
%     end
%     k = k+1;
% end
% Creconstructed_vec = nodeL;
% 
% while(size(Creconstructed_vec,1) ~= 1)
%     Creconstructed_vec = rowExpandMat(Creconstructed_vec',4)';
% end
% 
% Freconstructed_vec = (T'*Creconstructed_vec')'; %T^-1=T' since T is unitary
% Freconstructed = zeros(Nx,Ny);
% for j = 1:Ny
%     i = 1 + (j-1)*Nx;
% 
%     Freconstructed(:,j) = Freconstructed_vec(i:(i-1+Nx));
% end
% 
% % Plot
% %
% figure(21)
% ax1=subplot(1,2,1);
% plot( 1:(L-1),entanglementEntropies(1,:)./entanglementEntropies(2,:),'-*');
% title("Relative Entanglement Entropy versus lengthscale");
% xlabel("Bond l for lengthscale 2^{-l}");
% ylabel("S/S_{max}(l)");
% grid on;
% 
% ax2=subplot(1,2,2);
% plot( (1:(L-1)),(entanglementEntropies(1,:)),'-*');
% title("Entanglement Entropy at each bond");
% xlabel("Bond l for lengthscale 2^{-l}");
% ylabel("S(l)");
% grid on;
% 
% xticks(ax1, 1:(L-1));
% xticks(ax2, 1:(L-1));
% 
% figure(22)
% surf(Freconstructed-F)
% shading interp
% 
% figure(23)
% surf(Freconstructed)
% shading interp
% 
% std_error = sqrt(sum(sum((Freconstructed-F).^2))
% mean_absF = sum(sum(abs(F)))/(Nx*Ny)

 end