%% load data

prepSR;
data_s = load(fullfile(processed_lfp_dir, 'meanSub_2020-01-30_16-03-00.mat'));

%% do analysis

options = struct;
options.artifacts = [6074, 6077];

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

%% plot

plot_multitaper(mt_res);

%% after PCA analysis: plot/save PCA change

pca_change('mt_res.mat');