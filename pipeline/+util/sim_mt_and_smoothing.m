function [output_spectrograms, time, freqs] = sim_mt_and_smoothing(test_signal)
% Simulate the multitaper and smoothing used for LFP analysis.
% Theoretically we might be able to just use this function to do the analysis itself,
% but it was written to test the time resolution of the analysis using test signals.
%
% Input is a vector or an n_chans x n_time matrix of signals.
% Output is an n_freqs x n_windows x n_chans array of processed spectrograms.

Fs = 1000;  % sample rate

% make a struct for the data as if it were coming from a preprocessed recording
data_s = struct;
data_s.finalSampR = Fs;
if isvector(test_signal)
    test_signal = reshape(test_signal, 1, []);
end
n_chans = size(test_signal, 1);
data_s.meanSubFullTrace = test_signal;
data_s.info = struct;
data_s.info.noiseChannels = [];
data_s.info.Probe1Name = 'Probe1';
data_s.info.Probe1Channels = 1:n_chans;
data_s.info.Probe1Indicies = 1:n_chans;

% "high-res" multitaper parameters from the mt_* scripts
mt_opts = struct;
mt_opts.chans = struct('Probe1', 1:n_chans);
mt_opts.chan_names = strcat('C', string(1:n_chans));
mt_opts.window = 6;
mt_opts.padbase = 60;
mt_opts.winstep = 0.1;
mt_opts.save = false;

s = warning;
warning('off', 'ProektLab:unknownProbeWarn');
mt_res = multitaper_analysis(data_s, mt_opts);
warning(s);

% preprocessing steps
pp_opts = struct;
pp_opts.in_name = 'pxx';
pp_opts.name = 'pxx_rankord';
pp_opts.freq_sm_type = 'med';
pp_opts.freq_sm_span = 10;
pp_opts.time_sm_type = 'exp';
pp_opts.time_sm_span = 120; % seconds
pp_opts.norm_type = 'rank';

mt_res_pp = mt_preprocess(mt_res, pp_opts);

output_spectrograms = cell2mat(reshape(mt_res_pp.pxx_rankord, 1, 1, []));
time = mt_res_pp.time_grid;
freqs = mt_res_pp.freq_grid;

end

