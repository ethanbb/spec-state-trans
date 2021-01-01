function transitions = get_state_transitions(mt_res_files, nmf_res_file, dist_radius)
% Get a more detailed table of transitions as defined by discrete classes from an NMF analysis.
% Includes times in seconds, change distances based on NMF scores and distance percentiles
% based on comparing within each channel and across all channels.
% dist_radius is how far back and forward to look when computing distance across a transition, in s.

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

% get radius in samples
if ~exist('dist_radius', 'var') || isempty(dist_radius)
    dist_radius = 3;
end
mt_opts = mt_res.options;
radius_samps = dist_radius / mt_opts.winstep;

% get where each good segment starts and ends in post-multitaper samples
seg_ends = cumsum(seg_lengths);
seg_starts = [1, seg_ends(1:end-1) + 1];
segs = [seg_starts; seg_ends];
n_segs = length(segs);

% load classes and scores to identify transitions (just use 1st repetition)
if ischar(nmf_res_file)
    nmf_res_file = matfile(nmf_res_file);
end

classes = nmf_res_file.nmf_classes;
classes = classes{1};
scores = nmf_res_file.nmf_V;
scores = scores{1};

time_axis = nmf_res_file.time_axis;
time_segs = arrayfun(@(kS) time_axis(segs(1, kS):segs(2, kS)), 1:n_segs, 'uni', false);
chan_names = nmf_res_file.chan_names;
n_chans = length(chan_names);

chan_trans = cell(n_chans, 1);
for kC = 1:n_chans    
    trans_times = cell(n_segs, 1);
    dwell_times = cell(n_segs, 1);
    start_classes = cell(n_segs, 1);
    end_classes = cell(n_segs, 1);
    trans_dists = cell(n_segs, 1);
    
    for kS = 1:n_segs
        time = time_segs{kS}(:);
        chan_classes = classes{kC}(segs(1, kS):segs(2, kS));
        chan_scores = scores{kC}(segs(1, kS):segs(2, kS), :);
        
        b_trans = [false; diff(chan_classes) ~= 0];
        trans_inds = find(b_trans);
        n_trans = length(trans_inds);
        
        trans_times{kS} = time(b_trans);
        dwell_times{kS} = [trans_times{kS}(2:end); time(end)] - trans_times{kS};
        end_classes{kS} = chan_classes(trans_inds);
        start_classes{kS} = chan_classes(trans_inds - 1);
        
        start_scores = zeros(n_trans, size(chan_scores, 2));
        end_scores = zeros(n_trans, size(chan_scores, 2));
        for kT = 1:n_trans
            pre_inds = max(1, trans_inds(kT)-radius_samps):trans_inds(kT)-1;
            post_inds = trans_inds(kT):min(length(time), trans_inds(kT)+radius_samps-1);
            start_scores(kT, :) = mean(chan_scores(pre_inds, :));
            end_scores(kT, :) = mean(chan_scores(post_inds, :));
        end
        
        score_change = end_scores - start_scores;
        trans_dists{kS} = vecnorm(score_change, 2, 2);

        % try change in scores that are relevant to the transition
%         trans_inds = find(b_trans);
%         start_scores_pre = chan_scores(sub2ind(size(chan_scores), trans_inds-1, start_classes{kS}));
%         end_scores_pre = chan_scores(sub2ind(size(chan_scores), trans_inds-1, end_classes{kS}));
%         diff_pre = end_scores_pre - start_scores_pre;
%         
%         start_scores_post = chan_scores(sub2ind(size(chan_scores), trans_inds, start_classes{kS}));
%         end_scores_post = chan_scores(sub2ind(size(chan_scores), trans_inds, end_classes{kS}));
%         diff_post = end_scores_post - start_scores_post;
%         
%         diff_change = diff_post - diff_pre;
%         trans_dists{kS} = diff_change;
    end
    
    n_trans = sum(cellfun('length', trans_times));
    
    % make a table out of it
    chan_trans{kC} = table(...
        kC * ones(n_trans, 1), ...
        repmat(chan_names(kC), n_trans, 1), ...
        cell2mat(trans_times), ...
        cell2mat(dwell_times), ...
        cell2mat(trans_dists), ...
        cell2mat(start_classes), ...
        cell2mat(end_classes), ...
        'VariableNames', ...
        {'chan', 'chan_name', 'time', 'dwell_time', 'dist', 'start_class', 'end_class'});
    
    % get within-channel percentiles
    [~, order] = sort(chan_trans{kC}.dist);
    chan_trans{kC}.chan_pctile = zeros(n_trans, 1);
    chan_trans{kC}.chan_pctile(order) = (0.5:n_trans-0.5) / n_trans;
end

% combine all and get overall percentile
transitions = vertcat(chan_trans{:});
n_trans = height(transitions);
[~, order] = sort(transitions.dist);
transitions.pctile = zeros(n_trans, 1);
transitions.pctile(order) = (0.5:n_trans-0.5) / n_trans;

end
