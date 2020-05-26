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

% % chans based on 13-15-00 CSD:
% options.chans = [9, 39];
% options.chan_names = {'V1', 'MC'};

% 4 channels each in V1 and MC (spread out)
% artifact in chan 54 @ 270-272
options.chans = [6, 11, 22, 28, 38, 43, 55, 60];
options.chan_names = {'V1 super', 'V1 mid-super', 'V1 mid-deep', 'V1 deep', ...
                      'MC super', 'MC mid-super', 'MC mid-deep', 'MC deep'};

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

mt_res = multitaper_analysis(data_s, options);

%% preprocess

mt_preprocess(mt_res, pp_options);

%% plot

plot_options.filename = 'multitaper.fig';
plot_options.clim = [-4, 4];

plot_multitaper(mt_res, plot_options);