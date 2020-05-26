%% set dataset info

prepSR;

recdate = '2020-03-05';
time = '12-50-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

%% do low-resolution analysis first

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [];

% % chans based on 13-15-00 CSD:
% options.chans = [13, 58];
% options.chan_names = {'V1', 'MC'};

% 4 channels each in V1 and MC (spread out)
% artifact in chans 53-54 @ 5747
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

options.save = true;
options.savedir = savedir;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);


%% plot

plot_options = struct('savedir', savedir);

plot_multitaper(mt_res, plot_options);