%% Note - old analysis not currently in use. Current script is nnmf_clustering.m in parent directory.
%% Load processed multitaper results

prepSR;

recdate = '2020-02-06';
time = '16-01-00';
save_dir = fullfile(results_dir, recdate, time);

res_path = fullfile(save_dir, 'mt_res.mat');
res_mfile = matfile(res_path, 'Writable', true);

mt_opts = res_mfile.options; % necessary due to MatFile quirk
Fw = 1 / mt_opts.winstep; % "window rate"

chans = res_mfile.name;
chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false);
n_chans = numel(chans);

%% Do PCA on both channels together. Keep a relatively large # of components, say to 99% variance
% 
% pca_opts_both = struct;
% pca_opts_both.pxx_name = 'pxx_pp';
% pca_opts_both.name = 'pxx_pca_joint';
% pca_opts_both.thresh_type = 'cumvar';
% pca_opts_both.thresh = 99;
% [pc_data_both, fh] = mt_pca(res_path, pca_opts_both);

%% Remove NaNs and take subsample of frequencies (from all chans)

% nnmf_data = res_mfile.pxx_pca_joint;
data = res_mfile.pxx_pp;
freq_grid = res_mfile.freq_grid;

data_cat = horzcat(data{:});
b_nan = any(isnan(data_cat));
data_nonan = data_cat(:, ~b_nan);

time_ds_factor = 20;
freq_ds_factor = 10;

freqinds_ds = 1:freq_ds_factor:length(freq_grid);
timeinds_ds = 1:time_ds_factor:size(data_nonan, 2);
data_nonan_freqds = data_nonan(freqinds_ds, :);
data_nonan_ds = data_nonan_freqds(:, timeinds_ds);

data2use = data_nonan_ds;

%% Do NNMF on downsampled data

[U, V, p] = sp_nnmf(exp(data2use.'), 11, [], [], 500000);

%% Sort components according to peak frequency (since they should be sparse)

[~, colmax] = max(V);
[~, order] = sort(colmax);
U = U(:, order);
V = V(:, order);

%% Make distance matrix

% try not using all the components so distances aren't weird
%data2use = data2use(1:7, :);

[n_dims, n_pts] = size(data2use);
is_diag = logical(eye(n_pts));

% use 1/distance^2 as connection weights, so as to de-emphasize far connections

dist2s = squareform(pdist(data2use.', 'squaredeuclidean'));
conn_weights = 1 ./ dist2s;
conn_weights = conn_weights ./ max(conn_weights(~is_diag));

% What to do about diagonal? Don't want NNMF to prioritize matching a huge diagonal.
% Here I treat closeness b/w point and itself as the average closeness of 2 consecutive points
conn_weights(is_diag) = mean(diag(conn_weights, 1));

%conn_weights = exp(-dists.^2 / epsilon);

%% Do symmetric NNMF

U2 = symnnmf(conn_weights, 5);
