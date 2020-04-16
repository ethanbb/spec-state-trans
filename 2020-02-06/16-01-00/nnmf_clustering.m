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

pca_opts_both = struct;
pca_opts_both.pxx_name = 'pxx_pp';
pca_opts_both.name = 'pxx_pca_joint';
pca_opts_both.thresh_type = 'cumvar';
pca_opts_both.thresh = 99;
[pc_data_both, fh] = mt_pca(res_path, pca_opts_both);

%% Take subsample of points (from all chans) so that our N x N matrix doesn't blow up the computer

ds_factor = 20;

pca_data = res_mfile.pxx_pca_joint;
nnmf_inds = cellfun(@(d) 1:ds_factor:size(d, 2), pca_data, 'uni', false);
nnmf_data = cellfun(@(d, i) d(:, i), pca_data, nnmf_inds, 'uni', false);

nnmf_data_all = horzcat(nnmf_data{:});
nnmf_data_isvalid = all(~isnan(nnmf_data_all));
data2use = nnmf_data_all(:, nnmf_data_isvalid);

%% Make distance matrix

% start w/ z-score of pca data
data2use = normalize(data2use, 2);

% try not using all the components so distances aren't weird
%data2use = data2use(1:7, :);

[n_dims, n_pts] = size(data2use);

% use 1/distance as connection weights, so as to de-emphasize far connections

dists = squareform(pdist(data2use.'));
conn_weights = 1 ./ dists;
conn_weights(1:(n_pts+1):end) = 0; % avoid Infs on diagonal

%conn_weights = exp(-dists.^2 / epsilon);

%% Do symmetric NNMF

U = symnnmf(conn_weights, 5);