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

% chans based on 15-43-00 CSD:
options.chans = [9, 42];
options.chan_names = {'V1', 'MC'};

options.save = false;

mt_res_lores = multitaper_analysis(data_s, options);

%% preprocess

pp_options = struct;
pp_options.name = 'pxx_pp';
pp_options.freq_sm_type = 'med';
pp_options.freq_sm_span = 40;
pp_options.time_sm_type = 'exp';
pp_options.time_sm_span = 60;
pp_options.norm_type = 'log_z';

mt_res_lores = mt_preprocess(mt_res_lores, pp_options);

%% plot & save lo-res

plot_options = struct;
plot_options.pxx_name = 'pxx_pp';
plot_options.save = true;
plot_options.savedir = savedir;
plot_options.filename = 'multitaper_lores.fig';

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