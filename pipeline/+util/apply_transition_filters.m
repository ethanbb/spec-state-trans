function [tt_filtered, class_segs] = apply_transition_filters(tt_all, time_segs, class_segs, score_segs, filters)
% tt_all is a table of all transitions in one rat. For each segment of each channel
% individually, apply filters sequentially. See help for apply_filter_1seg.
% second output is a modified version of the cell of discrete classes, taking filtered-out
% transitions into account.

n_chans = length(score_segs);
tt_filtered = cell(n_chans, 1);

for kC = 1:n_chans
    tt_chan = tt_all(tt_all.chan == kC, :);
    chan_score_segs = score_segs{kC};
    chan_class_segs = class_segs{kC};
    
    n_segs = length(chan_score_segs);
    tt_filtered{kC} = cell(n_segs, 1);
    
    for kS = 1:n_segs
        tt_seg = tt_chan(tt_chan.segment == kS, :);
        seg_times = time_segs{kS};
        seg_scores = chan_score_segs{kS};
        seg_classes = chan_class_segs{kS};
        
        for kF = 1:length(filters)
            [tt_seg, seg_classes] = apply_filter_1seg(tt_seg, seg_times, seg_classes, seg_scores, filters(kF));
        end        
        tt_filtered{kC}{kS} = tt_seg;
        chan_class_segs{kS} = seg_classes;
    end
    
    tt_filtered{kC} = vertcat(tt_filtered{kC}{:});
    class_segs{kC} = chan_class_segs;
end
tt_filtered = vertcat(tt_filtered{:});

end

function [tt_seg, seg_classes] = apply_filter_1seg(tt_seg, seg_times, seg_classes, seg_scores, filter)
% given tt_seg, a table of transitions in one recording segment, along the the corresponding
% score matrix and timestamps in seconds, apply conditions in the struct filter to remove
% transitions while keeping dwell_time, start_class, and end_class consistent.
% assumes tt_seg is sorted in time.

curr_ind = 1;
while curr_ind <= height(tt_seg)
    curr_trans = tt_seg(curr_ind, :);
    have_next_trans = curr_ind ~= height(tt_seg);
    
    % keep track of whether we should consider the next transition as well
    can_merge_next = false;
    if have_next_trans
        next_trans = tt_seg(curr_ind+1, :);
        can_merge_next = true;
    end
    
    % get indices defining previous and current states
    state_samples = struct;
    if curr_ind > 1
        state_samples.prev_start = find(seg_times == tt_seg.time(curr_ind-1));
    else
        state_samples.prev_start = 1;
    end
    state_samples.curr_start = find(seg_times == curr_trans.time);
    state_samples.prev_end = state_samples.curr_start - 1;
    if have_next_trans
        state_samples.curr_end = find(seg_times == next_trans.time) - 1;
    else
        state_samples.curr_end = length(seg_times);
    end

    % if the transition is from a class to itself, always remove
    if curr_trans.start_class == curr_trans.end_class
        [tt_seg, curr_ind, seg_classes] = remove_transition(tt_seg, curr_ind, 'left', state_samples, seg_classes);
        continue;
    end
    
    % now check that all the conditions for removal are satisfied
    if ~isempty(filter.dist_pctile)
        if curr_trans.pctile >= filter.dist_pctile
            curr_ind = curr_ind + 1;
            continue;
        elseif have_next_trans
            can_merge_next = can_merge_next && next_trans.pctile < filter.dist_pctile;
        end
    end
    
    if ~isempty(filter.dist_chan_pctile)
        if curr_trans.chan_pctile >= filter.dist_chan_pctile
            curr_ind = curr_ind + 1;
            continue;
        elseif have_next_trans
            can_merge_next = can_merge_next && next_trans.chan_pctile < filter.dist_chan_pctile;
        end
    end
    
    % at this point current transition should be eliminated based on distance criteria.
    % check left and right state/segment separately.
    % if eligible, want to merge previous state to the "right".
    % otherwise, see whether we want to merge next state to left or right.
    can_merge_right = true;
    can_merge_left = true;
    tiebreak_doit = true;  % whether to do the merge left if there's a tie
    
    % check dwell time
    if ~isempty(filter.dwell_time)
        left_dwell = seg_times(state_samples.prev_end) - seg_times(state_samples.prev_start);
        can_merge_right = can_merge_right && left_dwell < filter.dwell_time;
        
        curr_dwell = seg_times(state_samples.curr_end) - seg_times(state_samples.curr_start);
        curr_dwell_below_thresh = curr_dwell < filter.dwell_time;
        can_merge_left = can_merge_left && curr_dwell_below_thresh;
        can_merge_next = can_merge_next && curr_dwell_below_thresh;
    end
    
    % check mean score ratio
    if ~isempty(filter.mean_score_ratio)
        prev_class = curr_trans.start_class;
        curr_class = curr_trans.end_class;
        
        scores_prev = seg_scores(state_samples.prev_start:state_samples.prev_end, :);
        mean_ratio_prev = mean(scores_prev(:, prev_class)) / mean(scores_prev(:, curr_class));
        can_merge_right = can_merge_right && mean_ratio_prev < filter.mean_score_ratio;
        
        scores_curr = seg_scores(state_samples.curr_start:state_samples.curr_end, :);
        mean_ratio_curr = mean(scores_curr(:, curr_class)) / mean(scores_curr(:, prev_class));
        can_merge_left = can_merge_left && mean_ratio_curr < filter.mean_score_ratio;
        
        if can_merge_next
            next_class = next_trans.end_class;
            mean_ratio_next = mean(scores_curr(:, curr_class)) / mean(scores_curr(:, next_class));
            can_merge_next = can_merge_next && mean_ratio_next < filter.mean_score_ratio;
        end
        
        if can_merge_next && strcmp(filter.tiebreaker, 'mean_score_ratio')
            tiebreak_doit = mean_ratio_curr < mean_ratio_next;
        end
    end
    
    % now merge if called for
    if can_merge_right
        [tt_seg, curr_ind, seg_classes] = remove_transition(tt_seg, curr_ind, 'right', state_samples, seg_classes);
        continue;
    end
    
    if can_merge_left && can_merge_next
        % resolve other tiebreakers
        switch filter.tiebreaker
            case 'dist_pctile'
                tiebreak_doit = curr_trans.pctile < next_trans.pctile;
            case 'dist_chan_pctile'
                tiebreak_doit = curr_trans.chan_pctile < next_trans.chan_pctile;
            case 'left'
                tiebreak_doit = true;
            case 'right'
                tiebreak_doit = false;
            case 'mean_score_ratio' % already handled
            otherwise
                error('Unknown tiebreaker: "%s"', filter.tiebreaker);
        end
        
        can_merge_left = tiebreak_doit;
    end
    
    if can_merge_left
        [tt_seg, curr_ind, seg_classes] = remove_transition(tt_seg, curr_ind, 'left', state_samples, seg_classes);
        continue;
    end
    
    curr_ind = curr_ind + 1;
end

end

function [tt_seg, next_ind, seg_classes] = remove_transition(tt_seg, ind, direction, state_samples, seg_classes)
% remove transition in position ind from the table. return the updated table and
% the index of the next transition that should be examined.
% direction should be 'left' or 'right' and specifies whether we take the class
% for the merged state from the left or right of the deleted transition.
% if the deletion creates a fake "transition" from a class to itself, it should be
% removed in a later call as apply_filter moves along the segment.

switch direction
    case 'left'
        new_class = tt_seg.start_class(ind);
        next_ind = ind;
        seg_classes(state_samples.curr_start:state_samples.curr_end) = new_class;
        mergeleft = true;
    case 'right'
        new_class = tt_seg.end_class(ind);
        next_ind = max(ind-1, 1);
        seg_classes(state_samples.prev_start:state_samples.prev_end) = new_class;
        mergeleft = false;
    otherwise
        error('Invalid merge "direction"');
end

if ind > 1
    tt_seg.dwell_time(ind-1) = sum(tt_seg.dwell_time([ind-1, ind]));
    if ~mergeleft
        tt_seg.end_class(ind-1) = new_class;
    end
end

if ind < height(tt_seg) && mergeleft
    tt_seg.start_class(ind+1) = new_class;
end

tt_seg(ind, :) = [];

end
