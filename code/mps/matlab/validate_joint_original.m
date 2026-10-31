% Original three-component error convention at the 64^3 grid.
[root,package]=mps_data_root();
addpath(fullfile(package,'original','code','PDE_matfiles','MPSBuilderNik'));
addpath(fullfile(package,'original','code','AnalysingResultsFromDNS','DNS','IncNs3DMatlab'));
target=fullfile(root,'recheck');
if ~exist(target,'dir'), mkdir(target); end
input=fullfile(target,'joint_original_inputs.mat');
if ~exist(input,'file'), input=fullfile(package,'reference','recheck','joint_original_inputs.mat'); end
load(input);
[chis,entropyU,entropyV,entropyW,spectrumU,spectrumV,spectrumW]=...
    findCritChi3D(U,V,W,[1e-4,1e-6,1e-8]);
output=fullfile(target,'joint_original_outputs.mat');
if exist('OCTAVE_VERSION','builtin')
    save('-mat7-binary',output,'chis','entropyU','entropyV','entropyW','spectrumU','spectrumV','spectrumW');
else
    save(output,'chis','entropyU','entropyV','entropyW','spectrumU','spectrumV','spectrumW','-v7');
end
disp(max(chis,[],1));
