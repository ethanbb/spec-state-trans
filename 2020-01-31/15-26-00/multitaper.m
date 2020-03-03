%% set dataset info

prepSR;

date = '2020-01-31';
time = '15-26-00';

savedir = fullfile(results_dir, date, time);

%% load and do analysis

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', date, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;

% try smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);


%% plot

plot_options = struct('savedir', savedir);

plot_multitaper(mt_res, plot_options);