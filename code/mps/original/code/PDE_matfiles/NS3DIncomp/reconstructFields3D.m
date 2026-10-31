%Reconstruct fluid fields

function [Ur, Vr, Wr] = reconstructFields3D(v1_TNT,v2_TNT, v3_TNT)

Ur = reconstructFromTNT(v1_TNT,8);
Vr = reconstructFromTNT(v2_TNT,8);
Wr = reconstructFromTNT(v3_TNT,8);

% Transpose back into Matlab convention: so x direction is along cols,
% while y is along rows.
% Ur = permute(Ur,[2,1,3]);
% Vr = permute(Vr,[2,1,3]);
% Wr = permute(Wr,[2,1,3]);

end