%% set dataset info

prepSR;

recdate = '2020-02-06';
time = '13-47-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

%% do low-resolution analysis first

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.artifacts = [
    1387, 1391
    5050, len_secs % (to end of recording)
    ];

% chans based on CSD:
options.chans = [9, 24, 41, 56];
options.chan_names = {'apical V1', 'basal V1', 'apical MC', 'basal MC'};

options.save = false;

mt_res_lores = multitaper_analysis(data_s, options);

%% plot & save lo-res

plot_options = struct;
plot_options.savedir = savedir;
plot_options.filename = 'multitaper_lores.fig';

plot_multitaper(mt_res_lores, plot_options);

%% do high-res analysis

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.artifacts = [
    1387, 1391
    5050, len_secs % (to end of recording)
    ];

% use different channels to match other recording on same day
options.chans = [18, 47];

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);


%% plot

plot_options = struct('savedir', savedir);

plot_multitaper(mt_res, plot_options);