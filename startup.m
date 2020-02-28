% To get paths to work on a new computer:
% - Set the path to include just the containing folder of this file (save
%   the pathdef.m file in the user folder, i.e. Documents/MATLAB)
% - Add all other folders of interest to the path by modifying this file.

% set some useful paths
global synology_dir
if ispc
    synology_dir = 'Z:';
else
    synology_dir = '/synology';
end

addpath(fullfile(synology_dir, 'Andi', 'Matlab'));
addpath(fullfile(synology_dir, 'brenna', 'eeglab2019_1'));
addpath(fullfile(synology_dir, 'brenna', 'States_rats', 'forEthan', 'misc'));