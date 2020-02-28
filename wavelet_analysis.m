function wavelet_analysis(data, chan, seg_bounds)
% Inputs:
%   data:       chans x samples matrix of LFP data
%   chan:       channel to analyze; must be <= size(data, 1)
%   seg_bounds: [start end] of segment of interest, in seconds (e.g. [0 1] is the first second)

Fs = 1000;
Fs_ds = 100;
ds_factor = Fs / Fs_ds;

% low-pass filter at 20 Hz
cutoff = 20;
b = fir1(Fs * 3 / cutoff, cutoff / (Fs / 2));
data_filt = filtfilt(b, 1, data(chan, :));
%data_filt = data;

% select segment of interest
seg_filt = data_filt(seg_bounds(1) * Fs + 1 : seg_bounds(2) * Fs);

% do wavelet transform
% wavelet paramters: [symmetry, time-bandwidth]. Decreasing time-bandwidth increases temporal precision.
[wt, f, coi] = cwt(seg_filt, Fs, ...
    'FrequencyLimits', [0.5, 20], ...
    'WaveletParamters', [3, 30], ...
    'VoicesPerOctave', 4);

% make scalogram
wt_abs = abs(wt);
wt_abs_norm = wt_abs ./ sum(wt_abs, 1); % normalize by total power
wt_log = log(wt_abs_norm);
wt_log_centered = wt_log - mean(wt_log, 2); % subtract mean normalized spectrum

% what we have now: "how is the distribution of power across frequencies (at each timestep)
% different from the average (in log space)?"

% put nans outside the cone of influence
wt_log_centered(f(:) < coi(:)') = nan;
%wt_abs(f(:) < coi(:)') = nan;

time = seg_bounds(1)+1/Fs:1/Fs:seg_bounds(2);

wt_log_centered_ds = downsample(wt_log_centered.', ds_factor).';
time_ds = downsample(time, ds_factor);

figure;
subplot(2, 1, 1);
surf(time_ds, f, wt_log_centered_ds, 'LineStyle', 'none');
%surf(time, f, wt_abs, 'LineStyle', 'none');
set(gca, 'YScale', 'log');
view(2);
colorbar;

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Log deviation from average frequency spectrum over time');

% amplitude velocity - let's use the non-normalized power
wt_abs_vel = diff(wt_abs, 1, 2) * Fs; % per second
wt_abs_speed = vecnorm(wt_abs_vel);
wt_abs_speed_ds = downsample(wt_abs_speed.', ds_factor).';
speed_time_ds = downsample(time(1:end-1) + 1/Fs/2, ds_factor);

subplot(2, 1, 2);
plot(speed_time_ds, wt_abs_speed_ds);
title('Speed of change in wavelet scalogram over time');
xlabel('Time (s)');

% phase? is there any interesting information in there?

end