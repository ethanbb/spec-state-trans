%% set dataset info

sr_dirs = prepSR;

recdate = '2020-03-05';
time = '12-50-00';

savedir = fullfile(sr_dirs.results, recdate, time);

data_s = load(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s_%s.mat', recdate, time)));

%% do low-resolution analysis first

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;

% artifact in chans 53-54 @ 5747

options.chans = [1, 13, 28, 33, 45, 60];
options.chan_names = {'V1_L23', 'V1_L4', 'V1_L5', 'M1_L23', 'M1_L4', 'M1_L5'};

options.save = false;

mt_res_lores = multitaper_analysis(data_s, options);

% plot to check
plot_options = struct;
plot_options.pxx_name = 'pxx';
plot_options.take_log = true;

plot_multitaper(mt_res_lores, plot_options);

%% do high-res analysis

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

options.save = true;
options.savedir = savedir;
options.filename = 'mt_res_layers.mat';

mt_res = multitaper_analysis(data_s, options);

%% plot

plot_multitaper(mt_res, plot_options);
