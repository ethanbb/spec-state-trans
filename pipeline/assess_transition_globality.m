function transitions = assess_transition_globality(transitions, radius)
% For each transition, find how many other transitions are within some radius (in s).
% For best results use a hard dwell time threshold when making the transition table.
% Radius should be <= dwell_threshold/2 to ensure at most 1 transition is found in each
% other channel near a particular transition.

trans_times = transitions.time;
chans = transitions.chan;
globality = arrayfun(@(tt, c) sum(chans ~= c & abs(trans_times - tt) < radius), trans_times, chans);
transitions.globality = globality + 1;
transitions.norm_globality = transitions.globality / n_chans / (2*radius);

end
