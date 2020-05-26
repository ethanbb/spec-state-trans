%% set dataset info

prepSR;

recdate = '2020-02-06';
time = '16-01-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

%% do low-resolution analysis first

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.artifacts = [
    832, 834
    ];

% chans 10-13 have artifact @ 224-226 s
% chans 12-14 have artifact @ 276-278 s
% avoid channels 15-17 due to artifact @ 694 seconds
% channel 14 is also bad - has long-persisting noise after artifact @ 832
% avoid channels 48-49 due to artifact @ 4533

% % chans based on 15-43-00 CSD:
% options.chans = [9, 42];
% options.chan_names = {'V1', 'MC'};

% 4 channels each in V1 and MC (spread out)
options.chans = [6, 9, 22, 28, 38, 43, 54, 60];
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

plot_multitaper(mt_res, plot_options);