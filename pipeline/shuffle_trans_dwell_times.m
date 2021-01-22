function [trans_shuffled, rng_state] = shuffle_trans_dwell_times(transitions)
% Given a table of transitions, for each segment of each channel, shuffle
% order of the transitions and make transition times reflect dwell times.

rng('shuffle');
rng_state = rng;

trans_shuffled = transitions;
n_chans = max(trans_shuffled.chan);
n_segs = max(trans_shuffled.segment);

for kC = 1:n_chans
    for kS = 1:n_segs
        trans_inds = find(trans_shuffled.chan == kC & trans_shuffled.segment == kS);
        if isempty(trans_inds)
            continue;
        end
        
        % Make sure it's sorted by time
        trans_shuffled(trans_inds, :) = sortrows(trans_shuffled(trans_inds, :), 'time');
        
        % First transition time stays the same; others are determined from it.
        start_time = trans_shuffled.time(trans_inds(1));
        dwell_times = trans_shuffled.dwell_time(trans_inds);
        shuffled_dwell_times = dwell_times(randperm(length(dwell_times)));
        trans_shuffled.dwell_time(trans_inds) = shuffled_dwell_times;
        trans_shuffled.time(trans_inds) = start_time + cumsum([0; shuffled_dwell_times(1:end-1)]);
    end
end

end
