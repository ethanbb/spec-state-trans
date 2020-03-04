%% set dataset info

prepSR;

date = '2020-02-06';
time = '16-01-00';

savedir = fullfile(results_dir, date, time);


%% load and do analysis

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', date, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.artifacts = [
    832, 834
    ];

% avoid channels 15-17 due to artifact @ 694 seconds
% channel 14 is also bad - has long-persisting noise after artifact @ 832
% avoid channels 48-49 due to artifact @ 4533
options.chans = [18, 47];

% try smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);


%% plot

plot_options = struct('savedir', savedir);

plot_multitaper(mt_res, plot_options);