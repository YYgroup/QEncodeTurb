function [critChis, entanEntrU, entanEntrV, SingValsU, SingValsV]  = findCritChi2D(U,V,desiredAccs)
%Finds the chi necessary to recast the 2D field (U,V) into MPS form at the accuracy
%desiredAcc.

%Prepare
N = length(U(1,:));
L = log2(N);
if(L ~= round(L))
    fprintf('\nError: make sure input F into findCritChi2D is a 2D function with both dimensions powers of 2. Exiting.\n');
    exit(1);
end
critChis = zeros(L-1,length(desiredAccs)); 

%Produce MPS decomposition
[SingValsU, entanEntrU,~] = createMPSandAnalyse2D(U,N);
[SingValsV, entanEntrV,~] = createMPSandAnalyse2D(V,N);

%Sequential search for appropriate chi (binary search would be overkill)
norm0 = sum(sum(U.^2+V.^2));
for l = 1:(L-1)
    
    for i=1:length(desiredAccs)
        error=1;
        chi_Guess = 1;
        desiredAcc=desiredAccs(i);
        while(error>desiredAcc)
            norm1 = sum((SingValsU(l,1:chi_Guess)).^2) + sum((SingValsV(l,1:chi_Guess)).^2);
            error = 1-(norm1/norm0);
            chi_Guess = chi_Guess+1;
        end
        %Here's the critical chi for the given accuracy
        critChis(l,i)=chi_Guess-1;
    end
end

end