%% set dataset info

prepSR;

recdate = '2020-03-11';
time = '14-32-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

%% do low-resolution analysis first

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    3039, 3041
    4120, 4169
    4194, 4211
    4222, 4225
    4736, 4738
    6109, 6111
];

options.chans = [8, 20, 32, 40, 52, 64];
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
