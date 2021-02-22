function transitions = assess_transition_globality(transitions, radius, n_chans_all)
% For each transition, find how many other transitions are within some radius (in s).
% For best results use a hard dwell time threshold when making the transition table.
% Radius should be <= dwell_threshold/2 to ensure at most 1 transition is found in each
% other channel near a particular transition.

if ~exist('n_chans_all', 'var')
    n_chans_all = length(unique(transitions.chan));
end

trans_times = transitions.time;
chans = transitions.chan;
globality = arrayfun(@(tt, c) sum(chans ~= c & abs(trans_times - tt) < radius), trans_times, chans);
transitions.globality = globality + 1;
transitions.norm_globality = transitions.globality / n_chans_all / (2*radius);

end
