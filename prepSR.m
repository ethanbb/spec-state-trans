global synology_dir;
global project_dir;
global raw_dir;
global processed_lfp_dir;

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