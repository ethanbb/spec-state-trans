function [transitions, mean_sync] = calc_transition_synchronization(transitions, n_chans_all)
% Given transition table for N channels, calculate the "SPIKE-Synchronization"
% measure for each transition, which is the fraction of other channels that have
% transitions close in time, defined adaptively based on the local inter-transition
% intervals (dwell times) of the channels involved.
% n_chans_all can be provided to override normalization in case some channels have no
% transitions at all.
% Source: https://journals.physiology.org/doi/full/10.1152/jn.00848.2014#_i8

transitions.sync_score = zeros(height(transitions), 1);

% treat each segment separately
segs = unique(transitions.segment);
if ~exist('n_chans_all', 'var')
    n_chans_all = length(unique(transitions.chan));
end

for kS = 1:length(segs)
    seg = segs(kS);
    b_seg = transitions.segment == seg;
    chans = unique(transitions.chan(b_seg));
    n_chans = length(chans);
    
    seg_trans_forchan = arrayfun(@(c) transitions(b_seg & transitions.chan == c, :), ...
        chans, 'uni', false);
    
    % find minimum nearby dwell times for each channel (to later calculuate tau)
    min_dwells_forchan = cell(n_chans, 1);
    for kC = 1:n_chans
        dwell_times = seg_trans_forchan{kC}.dwell_time(1:end-1); % only count ones between transitions
        min_dwells_forchan{kC} = min([inf; dwell_times], [dwell_times; inf]);
    end

    % now iterate over channels again and find the synchronization scores
    for kC = 1:n_chans
        % iterate through this channel's transitions, then other chans
        chan_trans = seg_trans_forchan{kC};
        n_this = height(chan_trans);
        sync_scores = zeros(n_this, 1);
        
        for kT = 1:n_this
            this_time = chan_trans.time(kT);
            
            for kO = 1:n_chans
                if kO == kC
                    continue;
                end
                
                % find nearest transition
                other_trans = seg_trans_forchan{kO};
                [nearest_dist, nearest_ind] = min(abs(other_trans.time - this_time));
                tau = min(min_dwells_forchan{kC}(kT), min_dwells_forchan{kO}(nearest_ind));
                
                if isfinite(tau) && nearest_dist < tau
                    sync_scores(kT) = sync_scores(kT) + 1/(n_chans_all-1);
                end               
            end
        end
        
        % put back into original table
        transitions.sync_score(b_seg & transitions.chan == chans(kC)) = sync_scores;
    end
end

mean_sync = mean(transitions.sync_score);
end
