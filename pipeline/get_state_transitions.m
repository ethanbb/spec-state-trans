function transitions = get_state_transitions(mt_res_files, nmf_res_file, options)
% Get a more detailed table of transitions as defined by discrete classes from an NMF analysis.
% Includes times in seconds, change distances based on NMF scores and distance percentiles
% based on comparing within each channel and across all channels.

opts = struct(...
    'dist_radius',  3,      ... seconds back and forward to look when computing dist across transition
    'filters',      struct( ... conjunctive conditions to remove transitions (remove if all are active)
                            ... this can be a struct array, in which case each element's
                            ... filters are applied sequentially. empty array = don't use condition.
            'dwell_time',       100,  ... merge adjoining segment A into B if A is < 20 seconds AND
            'mean_score_ratio', 1.1,  ... ratio of mean A vs. B class score in A segment is < 1.1 AND
            'dist_pctile',      [],   ... overall transition distance percentile is below this AND
            'dist_chan_pctile', [],   ... per-chan transition distance percentile is below this.
            'tiebreaker',       'dist_pctile' ... how to decide which way to merge (B -> A or B -> C)
                                              ... in case the conditions are satisfied either way.
        ) ...
    );

if exist('options', 'var') && isstruct(options)
    opts_in = fieldnames(options);
    for kO = 1:length(opts_in)
        opts.(opts_in{kO}) = options.(opts_in{kO});
    end
end

% in case of a single filter, allow specifying parameters directly
if isscalar(opts.filters)
    filter_params = fieldnames(opts.filters);
    for kP = 1:length(filter_params)
        if isfield(opts, filter_params{kP})
            opts.filters.(filter_params{kP}) = opts.(filter_params{kP});
        end
    end
end

% load seg windows
seg_lengths = cell(1, length(mt_res_files));
for kR = 1:length(mt_res_files)
    mt_res = mt_res_files{kR};
    if ischar(mt_res)
        mt_res = matfile(mt_res);
    end
    seg_lengths{kR} = cellfun('length', mt_res.seg_windows);
end
seg_lengths = cell2mat(seg_lengths);
n_segs = length(seg_lengths);

% get radius in samples
mt_opts = mt_res.options;
radius_samps = opts.dist_radius / mt_opts.winstep;

% get where each good segment starts and ends in post-multitaper samples
seg_ends = cumsum(seg_lengths);
seg_starts = [1, seg_ends(1:end-1) + 1];

% load classes and scores to identify transitions (just use 1st repetition)
if ischar(nmf_res_file)
    nmf_res_file = matfile(nmf_res_file);
end

classes = nmf_res_file.nmf_classes;
classes = classes{1};
scores = nmf_res_file.nmf_V;
scores = scores{1};

time_axis = nmf_res_file.time_axis;
time_axis = time_axis(:);

% divide things into "good data" segments
segmentize = @(var) arrayfun(@(ss, se) var(ss:se, :), seg_starts, seg_ends, 'uni', false);
time_segs = segmentize(time_axis);
class_segs = cellfun(segmentize, classes, 'uni', false);
score_segs = cellfun(segmentize, scores, 'uni', false);

chan_names = nmf_res_file.chan_names;
n_chans = length(chan_names);

all_trans = cell(n_chans, 1);
for kC = 1:n_chans    
    chan_trans = cell(n_segs, 1);
    
    for kS = 1:n_segs
        time = time_segs{kS};
        chan_classes = class_segs{kC}{kS};
        chan_scores = score_segs{kC}{kS};
        
        trans_inds = find(diff(chan_classes) ~= 0) + 1;
        n_trans = length(trans_inds);
        
        trans_times = time(trans_inds);
        dwell_times = [trans_times(2:end); time(end)] - trans_times;
        end_classes = chan_classes(trans_inds);
        start_classes = chan_classes(trans_inds - 1);
        
        start_scores = zeros(n_trans, size(chan_scores, 2));
        end_scores = zeros(n_trans, size(chan_scores, 2));
        for kT = 1:n_trans
            pre_inds = max(1, trans_inds(kT)-radius_samps):trans_inds(kT)-1;
            post_inds = trans_inds(kT):min(length(time), trans_inds(kT)+radius_samps-1);
            start_scores(kT, :) = mean(chan_scores(pre_inds, :));
            end_scores(kT, :) = mean(chan_scores(post_inds, :));
        end
        
        score_change = end_scores - start_scores;
        trans_dists = vecnorm(score_change, 2, 2);

        % make a table out of it
        chan_trans{kS} = table(...
            kC * ones(n_trans, 1), ...
            repmat(chan_names(kC), n_trans, 1), ...
            kS * ones(n_trans, 1), ...
            trans_times, dwell_times, ...
            trans_dists, start_classes, end_classes, ...
            zeros(n_trans, 1), zeros(n_trans, 1), ...
            'VariableNames', ...
            {'chan', 'chan_name', 'segment', 'time', 'dwell_time', ...
            'dist', 'start_class', 'end_class', 'chan_pctile', 'pctile'});
    end
    
    all_trans{kC} = vertcat(chan_trans{:});
    
    % get within-channel percentiles
    [~, order] = sort(all_trans{kC}.dist);
    n_trans = height(all_trans{kC});
    all_trans{kC}.chan_pctile(order) = (0.5:n_trans-0.5) / n_trans;
end

% combine all and get overall percentile
transitions = vertcat(all_trans{:});
[~, order] = sort(transitions.dist);
n_trans = height(transitions);
transitions.pctile(order) = (0.5:n_trans-0.5) / n_trans;

% finally, apply filters.
transitions = util.apply_transition_filters(transitions, time_segs, score_segs, opts.filters);

end

