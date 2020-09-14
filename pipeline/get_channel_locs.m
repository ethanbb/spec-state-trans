function channel_locs = get_channel_locs(data_s, chan_s, models)
% Get physical relative locations for selected channels on selected probes.
%
% Input data_s is a struct or matfile object obtained from loading a raw or
% preprocessed data file.
%
% Input chan_s is a struct of selected channels. See organize_lfp for accepted formats
% (the section regarding "when chans is a struct").
%
% If 'models' is provided, it should be a cell of the same length as the number of fields in chan_s,
% with each entry a string specifying the corresponding probe model, 
% or empty to to infer based on the probe name in the recording info (which is also
% the default if 'models' is not given).
%
% Output is struct with the same fields as chan_s. The values are k x ndims matrices,
% where k is the number of channels selected for the given probe and ndims
% is the dimensionality of the probe's channel map. Values are the relative position
% in each dimension in um, such that the distance in um between two channels i, j
% can be obtained with:
%
%   norm(channel_locs(i,:) - channel_locs(j,:))

ifo = data_s.info;

probe_names = fieldnames(chan_s);
n_probes = length(probe_names);

for kP = 1:n_probes
    name = probe_names{kP};
    req_chans = chan_s.(name);

    if exist('models', 'var') && ~isempty(models{kP})
        model = models{kP};
    else
        % annoying thing
        if startsWith(name, 'grid') || startsWith(name, 'ecog')
            model = ifo.ecogGridName;
        else
            model = ifo.([name, 'Name']);
        end
    end

    % get location map depending on probe model
    if startsWith(model, 'H3')
        % Cambridge Neurotech 64-channel
        loc_map = num2cell((1:64)' * 20);

    elseif startsWith(model, 'H4')
        % Cambridge Neurotech 32-channel
        loc_map = num2cell((1:32)' * 25);

    elseif startsWith(model, 'E64-500-20-60')
        % NeuroNexus ECoG array
        [xloc, yloc] = meshgrid((1:6) * 500, (1:11) * 500);

        % mark locations with no real electrode
        xloc(11, 3:4) = nan;
        yloc(11, 3:4) = nan;

        % combine
        loc_map = arrayfun(@(x, y) [y, x], xloc, yloc, 'uni', false);

    else
        error('Name %s does not match a known probe model', model);
    end

    % now get the entries corresponding to the selected channels
    if islogical(req_chans)
        req_locs = vertcat(loc_map{req_chans});
        % eliminate invalid channels (to match organize_lfp output)
        req_locs(any(isnan(req_locs), 2), :) = [];

    elseif isvector(req_chans)
        loc_map_valid = loc_map(cellfun(@(loc) all(~isnan(loc)), loc_map));
        req_locs = vertcat(loc_map_valid{req_chans});

    else
        % convert subscripts to linear indices first, then proceed as in logical indexing
        subs = num2cell(req_chans, 1);
        lin_inds = sub2ind(size(loc_map), subs{:});
        req_locs = vertcat(loc_map{lin_inds});
        % eliminate invalid channels
        req_locs(any(isnan(req_locs), 2), :) = [];
    end

    channel_locs.(name) = req_locs;
end

end
