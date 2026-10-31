% Execute the original constructor, analyzer and reconstructor unchanged.
[root,package]=mps_data_root();
addpath(fullfile(package,'original','code','PDE_matfiles','MPSBuilderNik'));
addpath(fullfile(package,'original','code','AnalysingResultsFromDNS','DNS','IncNs3DMatlab'));
addpath(fullfile(package,'original','code','AnalysingResultsFromDNS'));
target=fullfile(root,'validation');
if ~exist(target,'dir'), mkdir(target); end
input=fullfile(target,'original_inputs.mat');
if ~exist(input,'file'), input=fullfile(package,'reference','validation','original_inputs.mat'); end
load(input);
tic;
MPS_full=createMPS3D(F,64);
F_full=reconstructFromNik3D(MPS_full);
MPS_truncated=createMPS3D(F,16);
F_truncated=reconstructFromNik3D(MPS_truncated);
[original_spectra,original_entropy]=createMPSandAnalyse3D(F,64);
marker_mps=createMPS3D(marker,8);
marker_reconstructed=reconstructFromNik3D(marker_mps);
mem_fraction_full=memComp(8,64,5);
mem_fraction_truncated=memComp(8,16,5);
original_seconds=toc;
output=fullfile(target,'original_outputs.mat');
if exist('OCTAVE_VERSION','builtin')
    save('-mat7-binary',output,'MPS_full','F_full','MPS_truncated','F_truncated',...
        'original_spectra','original_entropy','marker_reconstructed',...
        'mem_fraction_full','mem_fraction_truncated','original_seconds');
else
    save(output,'MPS_full','F_full','MPS_truncated','F_truncated',...
        'original_spectra','original_entropy','marker_reconstructed',...
        'mem_fraction_full','mem_fraction_truncated','original_seconds','-v7');
end
fprintf('Original source validation completed in %.3f seconds.\n',original_seconds);
