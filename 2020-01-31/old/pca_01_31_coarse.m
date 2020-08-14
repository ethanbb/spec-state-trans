%% Do PCA on the 2 stim-free recordings on 01-31, after multitaper analysis with 6-second window

prepSR;

result_files = {
    fullfile(results_dir, '2020-01-31', '12-52-00', 'mt_res_coarse.mat')
    fullfile(results_dir, '2020-01-31', '15-26-00', 'mt_res_coarse.mat')
    };

% take components explaining at least 1% of variance
pca_opts = struct;
pca_opts.name = 'pxx_pca_01_31';
pca_opts.thresh_type = 'cumvar';
pca_opts.thresh = 75;

mt_pca(result_files, pca_opts);