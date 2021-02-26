function [transitions, mean_sync] = calc_transition_synchronization(transitions, n_chans_all)
% Given transition table for N channels, calculate the "SPIKE-Synchronization"
% measure for each transition, which is the fraction of other channels that have
% transitions close in time, defined adaptively based on the local inter-transition
% intervals (dwell times) of the channels involved.
% n_chans_all can be provided to override normalization in case some channels have no
% transitions at all.
% Source: https://journals.physiology.org/doi/full/10.1152/jn.00848.2014#_i8

transitions.sync_score = zeros(height(transitions), 1);
transitions.tau = zeros(height(transitions), 1);  % minimum of surrounding dwell times

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
    
    chan_trans_inds = arrayfun(@(c) find(b_seg & transitions.chan == c), chans, 'uni', false);
    
    % find minimum nearby dwell times for each channel (to later calculuate tau)
    for kC = 1:n_chans
        dwell_times = transitions.dwell_time(chan_trans_inds{kC}(1:end-1));
        transitions.tau(chan_trans_inds{kC}) = min([inf; dwell_times], [dwell_times; inf]);
    end

    % now iterate over channels again and find the synchronization scores
    for kC = 1:n_chans
        % iterate through this channel's transitions, then other chans
        n_this = length(chan_trans_inds{kC});
        
        for kT = 1:n_this
            this_ind = chan_trans_inds{kC}(kT);
            this_time = transitions.time(this_ind);
            
            for kO = 1:n_chans
                if kO == kC
                    continue;
                end
                
                % find nearest transition
                other_trans_times = transitions.time(chan_trans_inds{kO});
                [nearest_dist, nearest_ind] = min(abs(other_trans_times - this_time));
                tau = min(transitions.tau([this_ind, chan_trans_inds{kO}(nearest_ind)])) / 2;
                
                if isfinite(tau) && nearest_dist < tau
                    transitions.sync_score(this_ind) = transitions.sync_score(this_ind) + ...
                        1/(n_chans_all-1);
                end               
            end
        end
    end
end

mean_sync = mean(transitions.sync_score);
end
