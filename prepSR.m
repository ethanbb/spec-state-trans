global synology_dir;
global project_dir;
global raw_dir;
global processed_lfp_dir;
global code_dir;
global results_dir;

if isempty(synology_dir)
    if ispc
        synology_dir = 'Z:';
    else
        synology_dir = '/synology';
    end
end

project_dir = fullfile(synology_dir, 'brenna', 'States_rats');
raw_dir = fullfile(project_dir, 'RawData');
processed_lfp_dir = fullfile(project_dir, 'meanSubtracted_fullTrace');
code_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir, 'forEthan', 'spec-state-trans');

addpath(fullfile(synology_dir, 'Andi', 'Matlab'));
addpath(fullfile(synology_dir, 'brenna', 'eeglab2019_1'));
addpath(fullfile(code_dir, 'misc'));