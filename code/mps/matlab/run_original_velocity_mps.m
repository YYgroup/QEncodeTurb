function [MPS,relative_error] = run_original_velocity_mps(n,component,chi)
% Direct original-code reproduction, independent of all quantum modules.
% Example: run_original_velocity_mps(15,'ux',64)
% The authors' pointwise index loop is intentionally unchanged here.
[root,package]=mps_data_root();
addpath(fullfile(package,'original','code','PDE_matfiles','MPSBuilderNik'));
if ~ismember(n,[15,18,21,24,27]) || ~any(strcmp(component,{'ux','uy','uz'}))
    error('Invalid grid or component');
end
F=load_velocity_npy(fullfile(root,'fields',sprintf('n%d',n),[component,'.npy']));
MPS=createMPS3D(F,chi);
recovered=reconstructFromNik3D(MPS);
relative_error=norm(recovered(:)-F(:))/norm(F(:));
output=fullfile(root,'original_runs');
if ~exist(output,'dir'), mkdir(output); end
save(fullfile(output,sprintf('n%d_%s_chi%d.mat',n,component,chi)),
     'MPS','n','component','chi','relative_error','-v7');
fprintf('n=%d, %s, chi=%d, relative L2 error=%.12g\n',n,component,chi,relative_error);
end
