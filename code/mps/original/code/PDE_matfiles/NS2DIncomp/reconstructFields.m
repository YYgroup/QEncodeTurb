%Reconstruct fluid fields

function [Ur, Vr] = reconstructFields(v1_TNT,v2_TNT)

Ur = reconstructFromTNT(v1_TNT,4);
Vr = reconstructFromTNT(v2_TNT,4);

%L = log2(length(Ur));

%VisualiseFlowField(Ur,Vr);

end