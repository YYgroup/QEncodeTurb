function [data_root,package_root] = mps_data_root()
% Shared portable paths for the optional original MATLAB/Octave checks.
package_root=fileparts(fileparts(mfilename('fullpath')));
data_root=getenv('MPS_WORKDIR');
if isempty(data_root), data_root=fullfile(package_root,'generated'); end
if ~exist(data_root,'dir'), mkdir(data_root); end
end
