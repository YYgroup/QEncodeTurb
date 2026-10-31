%Reconstruct fields

function F = reconstructFromTNT(f_TNT,d)

f_nik = convertMPS_TNTNik(f_TNT);

if(d==2)
    F = reconstructFromNik1D(f_nik);
elseif(d == 4)
    F = reconstructFromNik2D(f_nik);
    % Transpose back into Matlab convention: so x direction is along cols,
    % while y is along rows.
    F = F';
elseif( d == 8)
        F = reconstructFromNik3D(f_nik);
        % Transpose back into Matlab convention: so x direction is along cols,
        % while y is along rows.
        F = permute(F,[2,1,3]);
else
    error("Local physical dimension must be either 4 (for 2D) or 8 (for 3D). Exiting.");
end

end