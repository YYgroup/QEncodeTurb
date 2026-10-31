%Converts a MPS from TNT form to Nik's form

function MPS_Nik = convertMPS_TNTNik(MPS_TNT)

L = length(MPS_TNT.nodes);

MPS_Nik = {};

%Every tensor in MPS_Nik must have its running indices ordered as D, L and 
%finally R in order of significance. Do this for every lengthscale.

for l = 1:L
    
    TNT_tensorVals = MPS_TNT.nodes(l).tensor.elems.vals;
    
    Dnum_TNT = MPS_TNT.nodes(l).indices{strfind(MPS_TNT.nodes(l).ids,'D')} + 1;
    Lnum_TNT = MPS_TNT.nodes(l).indices{strfind(MPS_TNT.nodes(l).ids,'L')} + 1;
    Rnum_TNT = MPS_TNT.nodes(l).indices{strfind(MPS_TNT.nodes(l).ids,'R')} + 1;
    
    TNT_dims =  MPS_TNT.nodes(l).tensor.dims;
    D_dim =  prod(TNT_dims(Dnum_TNT));
    L_dim =  prod(TNT_dims(Lnum_TNT));
    R_dim =  prod(TNT_dims(Rnum_TNT));
    
    %Need this because sometimes legs are "fat" and made up of multiple
    %"thin" ones, in TNT. Very annoying, but that's just how it is.   
    if(length(Dnum_TNT)>1)
        offset = length(Dnum_TNT)-1;
        Dnum_TNT=Dnum_TNT(1);
        if(Dnum_TNT<Lnum_TNT(1))
            Lnum_TNT = Lnum_TNT -offset;
        end
        if(Dnum_TNT<Rnum_TNT(1))
            Rnum_TNT = Rnum_TNT -offset;
        end
    end
    if(length(Lnum_TNT)>1)
        offset = length(Lnum_TNT)-1;
        Lnum_TNT=Lnum_TNT(1);
        if(Lnum_TNT<Dnum_TNT(1))
            Dnum_TNT = Dnum_TNT -offset;
        end
        if(Lnum_TNT<Rnum_TNT(1))
            Rnum_TNT = Rnum_TNT -offset;
        end
    end
    if(length(Rnum_TNT)>1)
        offset = length(Rnum_TNT)-1;
        Rnum_TNT=Rnum_TNT(1);
        if(Rnum_TNT<Dnum_TNT(1))
            Dnum_TNT = Dnum_TNT -offset;
        end
        if(Rnum_TNT<Lnum_TNT(1))
            Lnum_TNT = Lnum_TNT -offset;
        end
    end
    
    TNT_dims(Dnum_TNT) = D_dim;
    TNT_dims(Lnum_TNT) = L_dim;
    TNT_dims(Rnum_TNT) = R_dim;
    
    Nik_tensor = zeros(D_dim*L_dim,R_dim);
    
    for k = 1:R_dim
        for j = 1:L_dim
            for i = 1:D_dim
                index_Nik = (i-1) + D_dim*(j-1) + D_dim*L_dim*(k-1) +1; %Index in Niks tensor
                
                %indices_Nik = [i,j,k];
                indices_TNT = zeros(1,3);
                
                indices_TNT(Dnum_TNT) = i;
                indices_TNT(Lnum_TNT) = j;
                indices_TNT(Rnum_TNT) = k;         
                index_TNT = (indices_TNT(1)-1) + TNT_dims(1)*(indices_TNT(2)-1)...
                    + TNT_dims(1)*TNT_dims(2)*(indices_TNT(3)-1); %Corresponding index of the TNT tensor
                
                Nik_tensor(index_Nik) = TNT_tensorVals(index_TNT + 1);   
            end
        end
    end
    
    MPS_Nik{l} = Nik_tensor;
end



end
