function [critChis,entanEntrU,entanEntrV,entanEntrW,SingValsU,SingValsV,SingValsW]  = findCritChi3D(U,V,W,desiredAccs)
%Finds the chi necessary to recast the 2D field (U,V) into MPS form at the accuracy
%desiredAcc.

%Prepare
N = length(U(1,:,:));
L = log2(N);
if(L ~= round(L))
    fprintf('\nError: make sure input F into findCritChi3D is a 3D function with both dimensions powers of 2. Exiting.\n');
    exit(1);
end
critChis = zeros(L-1,length(desiredAccs)); 

%Produce MPS decomposition
[SingValsU,entanEntrU] = createMPSandAnalyse3D(U,ceil(sqrt(N^3)));
[SingValsV,entanEntrV] = createMPSandAnalyse3D(V,ceil(sqrt(N^3)));
[SingValsW,entanEntrW] = createMPSandAnalyse3D(W,ceil(sqrt(N^3)));

%Sequential search for appropriate chi (binary search would be overkill)
norm0 = sum(sum(sum(U.^2+V.^2 + W.^2)));
for l = 1:(L-1)
    
    for i=1:length(desiredAccs)
        error=1;
        chi_Guess = 1;
        desiredAcc=desiredAccs(i);
        while(error>desiredAcc)
            norm1 = sum( (SingValsU(l,1:chi_Guess)).^2) + sum((SingValsV(l,1:chi_Guess)).^2) + sum((SingValsW(l,1:chi_Guess)).^2);
            error = 1-(norm1/norm0);
            chi_Guess = chi_Guess+1;
        end
        %Here's the critical chi for the given accuracy
        critChis(l,i)=chi_Guess-1;
    end
end

end