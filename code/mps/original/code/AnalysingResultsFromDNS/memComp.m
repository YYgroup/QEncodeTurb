%Compares memory requirements of MPS representation vs exact. Inputs must
%be in the form (physical dimension, internal bond dimension, length
%of mps).

function memFrac = memComp(d,chi,L)
 
%% Number of parameters needed for exact representation

exact = d^L;

%% Calculate number of parameters needed for MPS-approximated representation

mpsApprox = 0;
L_dim = 1;
R_dim = 1;

for l = 1:L
   L_dim = R_dim;
   R_dim = min([d^l,d^(L-l), chi]);
   mpsApprox = mpsApprox + L_dim*R_dim*d;
   
   %Substract gauge DOFs
   if l<L
       mpsApprox = mpsApprox - R_dim^2;
   end
   
   %fprintf("\nat site %.0i, chi_L=%.0i and chi_R=%.0i",l,L_dim,R_dim);
end

memFrac = mpsApprox/exact;

end