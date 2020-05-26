%% set dataset info

prepSR;

recdate = '2020-03-11';
time = '12-31-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

%% do low-resolution analysis first

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    288, 290
];

% % chans based on 11-31-00 CSD (using 3 neighboring channels per region):
% options.chans = [20, 21, 22, 56, 57, 58];
% options.chan_names = {'V1_1', 'V1_2', 'V1_3', 'MC_1', 'MC_2', 'MC_3'};

% 4 channels each in V1 and MC (spread out)
options.chans = [6, 11, 22, 28, 38, 43, 54, 60];
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
plot_options.chans = [2, 5];

plot_multitaper(mt_res, plot_options);