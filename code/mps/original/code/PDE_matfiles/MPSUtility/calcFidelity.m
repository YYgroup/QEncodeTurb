%calculates the fidelity between two functions

function fidelity = calcFidelity(F,G)

if(length(F(1,1,:)) == 1) %assumes function is 2D
    
    overlap = sum(sum(F.*G));
    normF = sqrt(sum(sum(F.^2)));
    normG = sqrt(sum(sum(G.^2)));


else %otherwise 3D
   
    overlap = sum(sum(sum(F.*G)));
    normF = sqrt(sum(sum(sum(F.^2))));
    normG = sqrt(sum(sum(sum(G.^2))));
end

if(normF ==0 && normG == 0)
    fidelity = 1;
else
    fidelity = (overlap/(normF*normG))^2;
end

end