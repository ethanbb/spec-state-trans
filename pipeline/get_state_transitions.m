function transitions = get_state_transitions(mt_res_files, nmf_res_file)
% Get a more detailed table of transitions as defined by discrete classes from an NMF analysis.
% Includes times in seconds, change distances based on NMF scores and distance percentiles
% based on comparing within each channel and across all channels.

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
    start_classes = cell(n_segs, 1);
    end_classes = cell(n_segs, 1);
    trans_dists = cell(n_segs, 1);
    
    for kS = 1:n_segs
        time = time_segs{kS}(:);
        chan_classes = classes{kC}(segs(1, kS):segs(2, kS));
        chan_scores = scores{kC}(segs(1, kS):segs(2, kS), :);
        
        b_trans = [false; diff(chan_classes) ~= 0];
        trans_times{kS} = time(b_trans);
        end_classes{kS} = chan_classes(b_trans);
        start_classes{kS} = chan_classes(find(b_trans) - 1);
        
        end_scores = chan_scores(b_trans, :);
        start_scores = chan_scores(find(b_trans) - 1, :);
        score_change = end_scores - start_scores;
        trans_dists{kS} = vecnorm(score_change, 2, 2);
    end
    
    n_trans = sum(cellfun('length', trans_times));
    
    % make a table out of it
    chan_trans{kC} = table(...
        kC * ones(n_trans, 1), ...
        repmat(chan_names(kC), n_trans, 1), ...
        cell2mat(trans_times), ...
        cell2mat(trans_dists), ...
        cell2mat(start_classes), ...
        cell2mat(end_classes), ...
        'VariableNames', ...
        {'chan', 'chan_name', 'time', 'dist', 'start_class', 'end_class'});
    
    % get within-channel percentiles
    [~, order] = sort(chan_trans{kC}.dist);
    chan_trans{kC}.chan_pctile = (order - 0.5) / n_trans;
end

% combine all and get overall percentile
transitions = vertcat(chan_trans{:});
n_trans = height(transitions);
[~, order] = sort(transitions.dist);
transitions.pctile = (order - 0.5) / n_trans;

end
