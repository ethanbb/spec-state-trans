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

% % chans based on 16-33-00 CSD:
% options.chans = [18, 56];
% options.chan_names = {'V1', 'MC'};

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

plot_multitaper(mt_res, plot_options);