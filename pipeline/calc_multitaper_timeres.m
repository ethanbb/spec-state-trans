function time_res = calc_multitaper_timeres(Fs, win_length, step, n_tapers)
% Estimate the temporal resolution of multitaper by reference to the equivalent half-overlapping
% rectangular window with the same timestep, which we assume has temporal resolution = step size.
% The idea is that to detect a change in the signal (e.g. a step or the start of an oscillation),
% the maximum in the integrated amplitude envelope of the window from one step to the next must be
% large enough; otherwise more than one step is necessary. This means that the temporal resolution 
% will become higher as NW increases (while the frequency resolution decreases). Not sure if this is
% the right way to do it but there's some logic behind it.
%
% Fs: sample rate
% win_length, step: in seconds
% n_tapers: time-half-bandwidth product NW is inferred from this

NW = (n_tapers+1)/2;
win_samps = Fs * win_length;
step_samps = Fs * step;

[tapers, lambdas] = dpss(win_samps, NW);
tapers = tapers .* lambdas(:)';

% make equivalent rectangular window
rect_win_samps = 2*step_samps;
rect_win = ones(rect_win_samps, 1) / sqrt(rect_win_samps);

step_diffs_rect = get_onset_step_differences(rect_win, step_samps);
max_diff_rect = max(step_diffs_rect);

figure;
subplot(211);
plot(step_diffs_rect);
title('Rectangular window: diff in cumsum from 1 step');

% now figure out how many steps of the tapers we need to match it
max_diff_dpss = 0;
n_steps = 0;
while max_diff_dpss < max_diff_rect
    n_steps = n_steps + 1;
    step_diffs_dpss = get_onset_step_differences(tapers, step_samps * n_steps);
    max_diff_dpss = max(max(step_diffs_dpss));
end

time_res = step * n_steps;

subplot(212);
plot(step_diffs_dpss);
title(sprintf('Slepian tapers: diff in cumsum from %d steps', n_steps));

end

function step_diffs = get_onset_step_differences(windows, step_samps)
% compute the difference induced in the cumulative window sum from the onset of a signal occurring
% the given step_samps apart, at each point in the window.

if isvector(windows)
    windows = windows(:);
end
n_wins = size(windows, 2);

window_amp = abs(hilbert(windows));
window_amp_sum = cumsum(window_amp);
window_amp_sum_orig = [window_amp_sum; repmat(window_amp_sum(end, :), step_samps, 1)];
window_amp_sum_shifted = [zeros(step_samps, n_wins); window_amp_sum];

step_diffs = window_amp_sum_orig - window_amp_sum_shifted;

end
