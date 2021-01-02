function transitions = assess_transition_globality(transitions, pctile_threshold, dwell_threshold, radius)
% For each transition, find how many other transitions are within some radius (in s).
% First step ensures the min time between transitions in one channel is dwell_threshold (in s).
% Radius should be <= dwell_threshold/2 to ensure at most 1 transition is found in each
% other channel near a particular transition.

if radius > dwell_threshold / 2
    warning('Dwell threshold is too low - results may not be valid');
end

% Step 1: Filtering by distance percentile
transitions = transitions(transitions.pctile > pctile_threshold, :);

% Step 2: Condense transitions with low dwell time by taking only the one with highest distance
chans = unique(transitions.chan);
n_chans = length(chans);
chan_trans = cell(length(chans), 1);

for kC = 1:n_chans
    chan = chans(kC);
    chan_trans{kC} = sortrows(transitions(transitions.chan == chan, :), 'time');
    dwell_times = chan_trans{kC}.dwell_time;
    trans_below_thresh = find(dwell_times < dwell_threshold);
    
    inds_to_remove = [];
    curr_set = [];
    for kT = 1:length(trans_below_thresh)
        bt_ind = trans_below_thresh(kT);
        curr_set(end+1) = bt_ind; %#ok<AGROW>
        % find all consecutive transitions below the threshold before condensing
        if kT ~= length(trans_below_thresh) && trans_below_thresh(kT+1) == bt_ind + 1
            continue;
        end
        
        if bt_ind == height(chan_trans{kC})
            % just remove these ones, no following real transition
            inds_to_remove = [inds_to_remove, curr_set]; %#ok<AGROW>
        else
            % condense with following transition
            curr_set(end+1) = bt_ind + 1; %#ok<AGROW>
            % keep the one with maximum dist
            [~, trans_to_keep] = max(chan_trans{kC}.dist(curr_set));
            ind_to_keep = curr_set(trans_to_keep);
            inds_to_remove = [inds_to_remove, setdiff(curr_set, ind_to_keep)]; %#ok<AGROW>
            chan_trans{kC}.start_class(ind_to_keep) = chan_trans{kC}.start_class(curr_set(1));
            chan_trans{kC}.end_class(ind_to_keep) = chan_trans{kC}.end_class(curr_set(end));
            
            % fix dwell times
            pre_dwell = sum(chan_trans{kC}.dwell_time(curr_set(1:trans_to_keep-1)));
            post_dwell = sum(chan_trans{kC}.dwell_time(curr_set(trans_to_keep+1:end)));
            prev_real_trans = max(setdiff(1:curr_set(1)-1, inds_to_remove));
            if ~isempty(prev_real_trans)
                chan_trans{kC}.dwell_time(prev_real_trans) = chan_trans{kC}.dwell_time(prev_real_trans) + pre_dwell;
            end
            chan_trans{kC}.dwell_time(ind_to_keep) = chan_trans{kC}.dwell_time(ind_to_keep) + post_dwell;
            
            curr_set = [];
        end
    end
    chan_trans{kC}(inds_to_remove, :) = [];   
end

transitions = vertcat(chan_trans{:});

% Step 3: Find # of chans within radius for each transition
trans_times = transitions.time;
chans = transitions.chan;
globality = arrayfun(@(tt, c) sum(chans ~= c & abs(trans_times - tt) < radius), trans_times, chans);
transitions.globality = globality + 1;

end
