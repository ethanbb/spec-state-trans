function [V_shuffled, classes_shuffled, model, rseed] = shuffle_scores_markov(V, classes, transitions, model, shuffle)
% Shuffle NMF scores using discrete-time markov chain (for one bootstrap iteration)
% V is an N x k score matrix
% Classes is a vector of N discrete classes
% Transitions is a k x k matrix that tallies discrete class transitions
% Model input is optional - avoids recreating the markov chain each time if running in a loop.
% Shuffle - if false, don't shuffle the rng state.

if ~exist('shuffle', 'var') || shuffle
    rng('shuffle');
    rstate = rng;
    rseed = rstate.Seed;
end

if ~exist('model', 'var') || isempty(model)
    % Filter out classes that are never present
    valid_classes = find(sum(transitions, 2) > 0);
    trans_valid = transitions(valid_classes, valid_classes);
    
    model = dtmc(trans_valid, 'StateNames', valid_classes);
else
    valid_classes = str2double(model.StateNames);
end

% Initial state
start_state = 0;
% this loop should only repeat in the very edge case where the very last state
% in the vector is unique and is also randomly selected
while ~ismember(start_state, valid_classes)
    start_state = randsample(classes, 1);
end
x0 = zeros(1, length(valid_classes));
x0(valid_classes == start_state) = 1;

% Simulate
class_inds_shuffled = simulate(model, length(classes)-1, 'X0', x0);
classes_shuffled = valid_classes(class_inds_shuffled);

% Sample from real V to make shuffled V
V_shuffled = nan(size(V));
for k = valid_classes(:)'
    inds_byclass = find(classes == k);
    b_shuffled_inclass = classes_shuffled == k;
    V_shuffled_inds = randsample(inds_byclass, sum(b_shuffled_inclass), true);
    V_shuffled(b_shuffled_inclass, :) = V(V_shuffled_inds, :);
end
assert(all(all(~isnan(V_shuffled))), 'something''s wrong, not all timepoints assigned');

end
