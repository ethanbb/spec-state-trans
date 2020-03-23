function lfp_organized = organize_lfp(data_s)
% returns a nchan x ntime LFP matrix with channels in data_s reordered according to
% map in data_s.
%
% note that this switches probe 1 and probe 2 in order to put visual before
% motor, which is currently my convention.

% depending on whether we have raw or mean-subtracted/artifact-removed
% data, field containing data will be differed.

lfp_fields = {'meanSubFullTrace', 'meanSubData'}; % second one is for "snippits"
if any(isfield(data_s, lfp_fields))
    lfp_field = lfp_fields{isfield(data_s, lfp_fields)};
else
    warning('Input appears to be raw data - may contain artifacts, noise, etc.');
    lfp_field = 'LFPData';
    assert(isfield(data_s, lfp_field), 'LFP data not found in loaded dataset.');
end

lfp_organized = data_s.(lfp_field)([data_s.info.Probe2Indicies, data_s.info.Probe1Indicies], :, :);

end