function [lfp_organized, chans_used] = organize_lfp(data_s, req_chans, noise_inds, need_data, ...
                                                    b_snippits, time_inds)
% Returns a nchan x ntime LFP matrix with channels in data_s reordered according to
% map in data_s. The input can also be a matfile object.
%
% When req_chans is not passed or empty, returns all channels in order according to the "Channels"
% fields in the data file, which vary in the order they put the probes (but for the later ones,
% motor is before visual) and within each probe the order is superficial to deep, 'left' to 'right'.
%
% When req_chans is a double vector, this default order will be indexed by chans to produce the output
% channels.
%
% When req_chans is a struct, each field is a probe name (e.g. 'Probe1', 'Probe2', 'grid', 'fork') and
% the corresponding value is either:
%   * a vector of linear indices into the given probe
%   * a logical matrix the same shape as the channel map (e.g. 64x1 for H3 ('fork'), 32x1 for H4,
%     and 11x6 for 'grid', where the center 2 entries in the last row are ignored because they
%     do not correspond to real channels)
%   * a k by ndim matrix of indices into the channel map (like the inputs to sub2ind), in order to
%     get channels by subscript in a specific order. k is the number of channels requested, ndim
%     is the number of dimensions of the channel map.
%   * the string 'all' to get all channels on that probe (in default order according to channel map)
%
% If noise_inds is nonempty, it should be a vector of indices into the original
% (not reordered) data matrix that will be filtered out.
%
% If need_data is false, lfp_organized will be empty and only the second output will be returned.
%
% b_snippits is a boolean indicating whether 3D trial data is requested (default false).
%
% The first output is the concatenated channels selected from each probe, in the order
% the fields are listed in the struct (or the requested order if req_chans is a vector) 
%
% The second output is a struct, with each probe name as a field, containing which
% (reordered) channels were used. So if the 2 most superficial channels in a 32-channel linear probe
% were discarded, the corresponding struct entry will contain 3:32.
%
% Output will be the channels extracted from each entry, concatenated in order.

% depending on whether we have raw or mean-subtracted/artifact-removed
% data, field containing data will be different.

if isstruct(data_s)
    in_ds = @isfield;
    get_n_time = @(s, field, dim) size(s.field, dim);
else
    assert(isa(data_s, 'matlab.io.MatFile'), 'Invalid data_s');
    in_ds = @isprop;
    get_n_time = @size;
end

if ~exist('b_snippits', 'var') || isempty(b_snippits) || ~b_snippits
    b_snippits = false;
    if in_ds(data_s, 'meanSubFullTrace')
        lfp_field = 'meanSubFullTrace';
    else
        assert(in_ds(data_s, 'LFPData'), 'LFP data not found in loaded dataset.');
        warning('Input appears to be raw data - may contain artifacts, noise, etc.');
        lfp_field = 'LFPData';
    end
    n_time = get_n_time(data_s, lfp_field, 2);
else
    if in_ds(data_s, 'dataSnippits')
        lfp_field = 'dataSnippits';
    else
        assert(in_ds(data_s, 'meanSubData'), 'Snippit data not found in loaded dataset.');
        lfp_field = 'meanSubData';
    end
    n_time = get_n_time(data_s, lfp_field, 3);
end

if ~exist('time_inds', 'var') || isempty(time_inds)
    time_inds = 1:n_time;
end

% find what indices we have
ifo = data_s.info;
info_fields = fieldnames(ifo);
probe_names = regexp(info_fields, '^.*(?=Channels$)', 'match', 'once');
probe_names = probe_names(~cellfun('isempty', probe_names));
probe_names = setdiff(probe_names, {'noise'}); % since there's a field called noiseChannels
n_probes = length(probe_names);

% special case for grid because it's special
isgrid = startsWith(probe_names, 'grid') | startsWith(probe_names, 'ecog');

% make corresponding struct of absolute indices into lfp_field
probe_indices = struct; % indices into lfp_field
probe_indices_withinvalid = struct;
probe_channels = struct; % ordered channel number corresponding to each index

for kP = 1:n_probes
    name = probe_names{kP};
    if any(strcmp(info_fields, [name, 'Indicies']))
        indices = ifo.([name, 'Indicies']);
    elseif isgrid(kP)
        indices = ifo.gridIndicies;
    else
        error('Indices for %s not found', name);
    end
    
    valid_inds = indices ~= 0;
    lin_indices_valid = indices(valid_inds);  % remove non-existent channels
    
    chans = ifo.([name, 'Channels']);

    if ismember(1, lin_indices_valid)
        % indices are possibly relative to channels, rather than absolute
        lin_indices_valid = chans(lin_indices_valid);
        indices(valid_inds) = lin_indices_valid;
    end

    probe_indices.(name) = lin_indices_valid(:);
    probe_indices_withinvalid.(name) = indices;
    probe_channels.(name) = chans(:);
end
cat_indices = cell2mat(struct2cell(probe_indices));
cat_chans = cell2mat(struct2cell(probe_channels));
assert(length(unique(cat_indices)) == length(cat_indices), 'Probe channels overlap - something''s fishy');


% if no 'chans' provided, default to all
if ~exist('req_chans', 'var') || isempty(req_chans)
    req_chans = 1:length(cat_indices);
end

if isstruct(req_chans)
    req_probes = fieldnames(req_chans);
    n_req = length(req_probes);

    indices_to_use = cell(n_req, 1);
    for kP = 1:n_req
        probename = req_probes{kP};
        probe_chans = req_chans.(probename);

        % normalize probe name a bit
        if strncmp(probename, 'ecog', 4)
%             probename = 'grid';
            probename = 'ecog';
        end

        % choose indices based on type of req_chans
        if islogical(probe_chans)
            % index into precomputed indices map still in original shape (including invalid channels)
            probe_inds = probe_indices_withinvalid.(probename)(probe_chans);
            % eliminate any invalid channels
            probe_inds(probe_inds == 0) = [];

            indices_to_use{kP} = probe_inds(:);
            
        elseif ischar(probe_chans) || isstring(probe_chans)
            assert(strcmpi(probe_chans, 'all'), 'Unknown channel specifier %s', probe_chans);
            indices_to_use{kP} = probe_indices.(probename);
            
        elseif isvector(probe_chans)
            indices_to_use{kP} = probe_indices.(probename)(probe_chans);
        else
            % convert subscripts to linear indices first, then proceed as in logical indexing
            subs = num2cell(probe_chans, 1);
            lin_inds = sub2ind(size(probe_indices_withinvalid.(probename)), subs{:});

            probe_inds = probe_indices_withinvalid.(probename)(lin_inds);
            % eliminate invalid channels
            probe_inds(probe_inds == 0) = [];

            indices_to_use{kP} = probe_inds(:);
        end
    end
    indices_to_use = cell2mat(indices_to_use);

else
    % grab channels according to the ordered channel numbers
    
    % we have a vector of the index of each channel.
    % need to reorder so the channels are in increasing order, then index into it
    chan_inds_ordered(cat_chans) = cat_indices;
    indices_to_use = chan_inds_ordered(req_chans);
end

% take noise_inds into account
if exist('noise_inds', 'var') && ~isempty(noise_inds)
    indices_to_use = setdiff(indices_to_use, noise_inds, 'stable');
end

% get channels that we actually have left
chans_used = struct;
for kP = 1:n_probes
    name = probe_names{kP};
    this_inds = probe_indices.(name);
    chans = cell2mat(arrayfun(@(i) find(i == this_inds), indices_to_use(:), 'uni', false));
    chans_used.(name) = chans;
end

if ~exist('need_data', 'var') || isempty(need_data) || need_data
    if isscalar(indices_to_use)
        if b_snippits
            lfp_organized = data_s.(lfp_field)(indices_to_use, :, time_inds);
        else
            lfp_organized = data_s.(lfp_field)(indices_to_use, time_inds);
        end
    else
        lfp_all = data_s.(lfp_field); % for mfile
        if b_snippits
            lfp_organized = lfp_all(indices_to_use, :, time_inds);
        else
            lfp_organized = lfp_all(indices_to_use, time_inds);
        end
    end
else
    lfp_organized = [];
end

end
