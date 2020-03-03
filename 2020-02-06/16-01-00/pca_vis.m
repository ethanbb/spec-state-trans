%% PCA just on this recording, just for the visual cortex channel

prepSR;

date = '2020-02-06';
time = '16-01-00';
res_file = fullfile(results_dir, date, time, 'mt_res.mat');

pca_opts = struct;
pca_opts.name = 'pxx_pca_1rec_vis';
pca_opts.chans = 1;
pca_opts.thresh_type = 'comps';
pca_opts.thresh = 10;

mt_pca(res_file, pca_opts);