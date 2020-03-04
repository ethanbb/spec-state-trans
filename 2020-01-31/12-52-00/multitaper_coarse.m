%% set dataset info

prepSR;

date = '2020-01-31';
time = '12-52-00';

savedir = fullfile(results_dir, date, time);

%% load and do analysis

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', date, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.filename = 'mt_res_coarse.mat';

mt_res = multitaper_analysis(data_s, options);


%% plot

plot_options = struct;
plot_options.savedir = savedir;
plot_options.filename = 'multitaper_coarse.fig';

plot_multitaper(mt_res, plot_options);