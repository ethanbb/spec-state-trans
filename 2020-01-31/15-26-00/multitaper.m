%% load data

prepSR;
data_s = load(fullfile(processed_lfp_dir, 'meanSub_2020-01-31_15-26-00.mat'));

%% do analysis

options = struct;

% smaller window
options.window = 6;
options.winstep = 0.1;
options.padbase = 60;

mt_res = multitaper_analysis(data_s, options);

%% plot

plot_multitaper(mt_res);