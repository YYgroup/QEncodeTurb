%Converts a MPS from Nik's form to TNT form

function wf = convertMPS_NikTNT(MPS_Nik,d)

L = length(MPS_Nik);

%d = 8 in 3D, 4 in 2D.

for loop = 1:L
    %Feed the data into the TNT MPS
    tensor.elems_type = 'values';
        
    tensor.elems.vals = MPS_Nik{loop};
    tensor.qn_info.qn_dir = [0,0,0];
    tensor.qn_info.qn_index = {[],[],[]};
    
    %Now properly define the legs
    [chi_L,chi_R] = size(MPS_Nik{loop}); chi_L = chi_L/d; 
    
    tensor.dims = [d,chi_L, chi_R];
    
    wf.nodes(loop).tensor = tensor;
    wf.nodes(loop).ids = 'DLR'; %LRD is default
    wf.nodes(loop).indices = {0,1,2};
end

%Define starting and ending nodes
wf.start = 0;
wf.start_leg = 'L';
wf.end = L-1;
wf.end_leg = 'R';

%Build up connections
wf.connections = cell(L,L);
for loop=2:L
    wf.connections{loop-1,loop} = 'RL';
end

end
