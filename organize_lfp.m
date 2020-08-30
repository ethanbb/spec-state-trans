function lfp_organized = organize_lfp(data_s, chans)
% Returns a nchan x ntime LFP matrix with channels in data_s reordered according to
% map in data_s. The input can also be a matfile object.
%
% Note that this switches probe 1 and probe 2 in order to put visual before
% motor, which is currently my convention.
%
% When 'chans' is not passed or empty, returns all channels in default order, which is Probe2 (V1)
% before Probe1 (M1) or 'grid' before 'fork', and within each probe channels are in matrix order
% according to the channel map which is superficial to deep, 'left' to 'right'.
%
% When chans is a double vector, this default order will be indexed by chans to produce the output
% channels.
%
% When chans is a struct, each field is a probe name (e.g. 'Probe1', 'Probe2', 'grid', 'fork') and
% the corresponding value is either:
%   * a vector of linear indices into the given probe
%   * a logical matrix the same shape as the channel map (e.g. 64x1 for H3 ('fork'), 32x1 for H4,
%     and 11x6 for 'grid', where the center 2 entries in the last row are ignored because they
%     do not correspond to real channels)
%   * a k by ndim matrix of indices into the channel map (like the inputs to sub2ind), in order to
%     get channels by subscript in a specific order. k is the number of channels requested, ndim
%     is the number of dimensions of the channel map.
% The output is the concatenated channels selected from each probe, in the order the fields are
% listed in the struct.
%
% Output will be the channels extracted from each entry, concatenated in order.

% depending on whether we have raw or mean-subtracted/artifact-removed
% data, field containing data will be different.

if isstruct(data_s)
    in_ds = @isfield;
else
    assert(isa(data_s, 'matlab.io.MatFile'), 'Invalid data_s');
    in_ds = @isprop;
end

if in_ds(data_s, 'meanSubFullTrace')
    lfp_field = 'meanSubFullTrace';
elseif in_ds(data_s, 'meanSubData') % for "snippits"
    lfp_field = 'meanSubData';
else
    assert(in_ds(data_s, 'LFPData'), 'LFP data not found in loaded dataset.');
    warning('Input appears to be raw data - may contain artifacts, noise, etc.');
    lfp_field = 'LFPData';
end

% find what indices we have
ifo = data_s.info;
info_fields = fieldnames(ifo);
probe_names = regexp(info_fields, '^.*(?=Indicies$)', 'match', 'once');
probe_names = probe_names(~cellfun('isempty', probe_names));
n_probes = length(probe_names);

% make corresponding struct of absolute indices into lfp_field
probe_indices = struct;
probe_indices_withinvalid = struct;
for kP = 1:n_probes
    name = probe_names{kP};
    indices = ifo.([name, 'Indicies']);
    valid_inds = indices ~= 0;
    lin_indices_valid = indices(valid_inds);  % remove non-existent channels

    if ismember(1, lin_indices_valid)
        % indices are possible relative to channels, rather than absolute
        lin_indices_valid = ifo.([strrep(name, 'grid', 'ecog'), 'Channels'])(lin_indices_valid);
        indices(valid_inds) = lin_indices_valid;
    end

    probe_indices.(name) = lin_indices_valid(:);
    probe_indices_withinvalid.(name) = indices;
end
cat_indices = cell2mat(struct2cell(probe_indices));
assert(length(unique(cat_indices)) == length(cat_indices), 'Probe channels overlap - something''s fishy');


% if no 'chans' provided, default to all
if ~exist('chans', 'var') || isempty(chans)
    chans = 1:length(cat_indices);
end

if isstruct(chans)
    req_probes = fieldnames(chans);
    n_req = length(req_probes);

    indices_to_use = cell(n_req, 1);
    for kP = 1:n_req
        probename = req_probes{kP};
        req_chans = chans.(probename);

        % normalize probe name a bit
        if strncmp(probename, 'ecog', 4)
            probename = 'grid';
        end

        % choose indices based on type of req_chans
        if islogical(req_chans)
            % index into precomputed indices map still in original shape (including invalid channels)
            req_inds = probe_indices_withinvalid.(probename)(req_chans);
            assert(all(req_inds > 0), 'Channel selection for %s included invalid channels', probename);

            indices_to_use{kP} = req_inds(:);

        elseif isvector(req_chans)
            indices_to_use{kP} = probe_indices.(probename)(req_chans);
        else
            % convert subscripts to linear indices first, then proceed as in logical indexing
            subs = num2cell(req_chans, 1);
            lin_inds = sub2ind(size(probe_indices_withinvalid.(probename)), subs{:});

            req_inds = probe_indices_withinvalid.(probename)(lin_inds);
            assert(all(req_inds > 0), 'Channel selection for %s included invalid channels', probename);

            indices_to_use{kP} = req_inds(:);
        end
    end
    indices_to_use = cell2mat(indices_to_use);

else
    % this is all legacy crap that doesn't make much sense. use the struct form if possible.

    % put them in the default order
    probe2pos = find(strcmp(probe_names, 'Probe2'), 1); % for V1
    if ~isempty(probe2pos)
        probe_names = circshift(probe_names, 1 - probe2pos);
    end

    gridpos = find(strcmp(probe_names, 'grid'), 1);
    if ~isempty(gridpos)
        probe_names = circshift(probe_names, 1 - gridpos);
    end

    cat_indices = cell(n_probes, 1);
    for kP = 1:n_probes
        cat_indices{kP} = probe_indices.(probe_names{kP});
    end
    cat_indices = cell2mat(cat_indices);

    indices_to_use = cat_indices(chans);
end

lfp_all = data_s.(lfp_field); % for mfile
lfp_organized = lfp_all(indices_to_use, :, :);

end
