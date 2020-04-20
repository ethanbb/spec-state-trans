%% Load processed multitaper results

prepSR;

redo_all = false; % change to true to redo everything

redo_pp = false; % change to true to redo preprocessing

recdate = '2020-02-06';
time = '13-47-00';

res_path = fullfile(results_dir, recdate, time, 'mt_res.mat');

res_mfile = matfile(res_path, 'Writable', true);

mt_opts = res_mfile.options; % necessary due to MatFile quirk
Fw = 1 / mt_opts.winstep; % "window rate"

chans = res_mfile.name;
chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false);
n_chans = numel(chans);

%% Re-preprocess data with parameters we want and concatenate

freq_ds_factor = 10;

if ~isprop(res_mfile, 'pxx_pp') || redo_pp || redo_all
    pp_options = struct;
    pp_options.in_name = 'pxx';
    pp_options.name = 'pxx_pp';
    pp_options.freq_sm_type = 'med';
    pp_options.freq_sm_span = freq_ds_factor;
    pp_options.time_sm_type = 'exp';
    pp_options.time_sm_span = 60;
    pp_options.norm_type = 'log_z';
    
    all_data = mt_preprocess(res_mfile, pp_options);
else
    all_data = res_mfile.pxx_pp;
end

% concatenate and collect information needed to break back apart
b_nan = cell(n_chans, 1);
chan_inds = all_data;
chan_inds_nonan = all_data;

last_ind = 0;
last_ind_nonan = 0;
    
for kC = 1:n_chans
    chan_inds{kC} = last_ind + (1:size(all_data{kC}, 2));
    last_ind = chan_inds{kC}(end);

    b_nan{kC} = any(isnan(all_data{kC}));
    chan_inds_nonan{kC} = last_ind_nonan + (1:size(all_data{kC}, 2)-sum(b_nan{kC}));
    last_ind_nonan = chan_inds_nonan{kC}(end);
end

data_cat = horzcat(all_data{:}).';
b_nan_cat = horzcat(b_nan{:}).';

%% Remove NaNs and take subsample of frequencies (from all chans)

freq_grid = res_mfile.freq_grid;
freqinds_ds = 1:freq_ds_factor:length(freq_grid);
freq_grid_ds = freq_grid(freqinds_ds);

data_nonan_freqds = data_cat(~b_nan_cat, freqinds_ds);
data2use = data_nonan_freqds;

%% Clip data below 0.1% and above 99.9% quantiles to deal with extreme outliers
qs = quantile(data2use(:), [0.001, 0.999]);
data2use(data2use < qs(1)) = qs(1);
data2use(data2use > qs(2)) = qs(2);

%% Test how many NNMF components satisfactorily reduce error - on downsampled data

n_comp_options = 1:18;

data2use_ds = data2use(1:20:end, :);
BCV_err = nnmf_k_ver(exp(data2use_ds), round(size(data2use_ds) / 10), 10, ...
    n_comp_options(1), n_comp_options(end), [], 500000);

figure;
scatter(reshape(repmat(n_comp_options, size(BCV_err, 1), 1), [], 1), BCV_err(:), 'k', 'filled');
hold on;
plot(n_comp_options, mean(BCV_err), 'b', 'LineWidth', 1);
xticks(n_comp_options);
xlabel('Number of NMF components');
ylabel('Reconstruction error');
title('NMF parameter selection using cross-validation');

%% Do NNMF

% based on previous analysis, # of comps selected empirically based on where
% error seems to enter more of a linearly decreasing regime.
n_comps = 7;
[U_combined, V_combined, p] = sp_nnmf(exp(data2use), n_comps, [], [], 500000);

%% Sort components according to peak frequency (since they should be sparse) and plot

[~, colmax] = max(V_combined);
[~, order] = sort(colmax);
U_combined = U_combined(:, order);
V_combined = V_combined(:, order);

figure;
sanePColor(1:n_comps, freq_grid_ds, V_combined, false, true);
set(gca, 'YScale', 'log');
xticks(1:n_comps);
xlabel('Component #');
ylabel('Frequency (Hz)');
title('NMF components');

%% Split result back up into individual datasets

U_by_chan = cell(size(all_data));

for kC = 1:n_chans
    U_by_chan{kC} = nan(length(b_nan{kC}), n_comps);
    U_by_chan{kC}(~b_nan{kC}, :) = U_combined(chan_inds_nonan{kC}, :);
end

%% Plot component-space representations

% don't have to plot all of them... specify kP and kC then run this section.
figure;
sanePColor(res_mfiles{kP}.time_grid, 1:n_comps, U_by_chan{kC, kP}');
xlabel('Time (s)');
ylabel('Component #');
title(sprintf('NMF component representation of %s on %s', chans{kC}, recs{kP}));

%% Plot reconstructions

% again, don't have to plot all of them
figure;
sanePColor(res_mfiles{kP}.time_grid, freq_grid_ds, log(U_by_chan{kC, kP} * V_combined')', false, true);
set(gca, 'YScale', 'log');
set(gca, 'CLim', [-4, 4]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('NMF reconstruction of %s on %s (log space)', chans{kC}, recs{kP}));

%% Fit a GMM to component representation

% plausible #s of clusters - at least 4, maybe as many as ~10? Could be higher if several are
% qualitatively similar (may be multiple "clusters" for a single "class")
ncs = 4:10;
n_ncs = length(ncs);

models = cell(length(n_ncs), 1);
for kK = 1:n_ncs
    opts = statset('MaxIter', 400);
    % add regularization value to try to avoid ill-conditioned covariance matrices
    models{kK} = fitgmdist(U_combined, ncs(kK), 'RegularizationValue', 0.01, 'Replicates', 3, ...
        'Options', opts);
end

AICs = cellfun(@(m) m.AIC, models);
figure;
plot(ncs, AICs);
xlabel('# of clusters');
ylabel('AIC value');
title('GMM model selection (Akaike information criterion)');

%% Looks like 6 clusters may be a good tradeoff - take a look at it

best_ncs = 6;

best_model = models{ncs == best_ncs};

%% Do initial breakdown of clusters by channels, for next step

cluster_probs_temp = posterior(best_model, U_combined);
[~, assigned_clusts_temp] = max(cluster_probs_temp, [], 2);

clusts_by_chan_temp = cell(size(all_data));


time_grid = res_mfile.time_grid;

for kC = 1:n_chans
    clusts_by_chan_temp{kC} = nan(length(b_nan{kC}), 1);
    clusts_by_chan_temp{kC}(~b_nan{kC}, :) = assigned_clusts_temp(chan_inds_nonan{kC}, :);
end

%% Make transition frequency and times matrices
transmat = zeros(best_ncs);
    
for kC = 1:n_chans
    clusts = clusts_by_chan_temp{kC};
    for kT = 1:length(clusts)-1
        if ~isnan(clusts(kT)) && ~isnan(clusts(kT+1))
           transmat(clusts(kT), clusts(kT+1)) = transmat(clusts(kT), clusts(kT+1)) + 1;
        end
    end
end

% exclude self-loops to focus on transitions
transmat(1:best_ncs+1:end) = 0;

clear *_temp;

%% Sort clusters based on frequency of transition between them, using hierarchical clustering

% let proximity = max transition probability in either direction - don't want to discount
% one-way connections. and then distance = 1 / proximity, with eps to prevent infs
sym_prox = max(transmat, transmat.');
pseudo_dist = 1 ./ (squareform(sym_prox) + eps);

% hierarchically cluster and get permutation
z = linkage(pseudo_dist);
clust_perm = optimalleaforder(z, pseudo_dist);

% permute everything
best_model = gmdistribution(best_model.mu(clust_perm, :), best_model.Sigma(:, :, clust_perm), ...
    best_model.ComponentProportion(clust_perm));
transmat = transmat(clust_perm, clust_perm);

%% See what the centroids look like in NMF and frequency space
figure;
sanePColor(1:best_ncs, 1:n_comps, best_model.mu');
xlabel('GMM cluster');
ylabel('NMF component');
title('Component-cluster mapping');

mu_freqspace = V_combined * best_model.mu';
figure;
sanePColor(1:best_ncs, freq_grid_ds, mu_freqspace, false, true);
set(gca, 'YScale', 'log');
xlabel('GMM cluster');
ylabel('Frequency (Hz)');
title('Frequency-cluster mappping');

%% Plot a Markov model of class/state transitions

mc = dtmc(transmat);
figure;
graphplot(mc, 'ColorEdges', true);

%% Get most likely clusters and break into individual channels

cluster_probs = posterior(best_model, U_combined);
[~, assigned_clusts] = max(cluster_probs, [], 2);

clusts_by_chan = cell(size(all_data));
all_transitions = cell([length(all_data), size(transmat)]);

for kC = 1:n_chans
    clusts_by_chan{kC} = nan(length(b_nan{kC}), 1);
    clusts_by_chan{kC}(~b_nan{kC}, :) = assigned_clusts(chan_inds_nonan{kC}, :);

    clusts = clusts_by_chan{kC};
    for kT = 1:length(clusts)-1
        if ~isnan(clusts(kT)) && ~isnan(clusts(kT+1))
           all_transitions{kC, clusts(kT), clusts(kT+1)}(end+1) = time_grid(kT+1);
        end
    end
end

%% see the distances betwen centroids
cent_dist = squareform(pdist(best_model.mu));
figure;
sanePColor(1:best_ncs, 1:best_ncs, cent_dist);
set(gca, 'YDir', 'reverse');
title('Pairwise Euclidean distances between GMM centroids');
colorbar;

%% save transitions to file for further analysis

res_mfile.freq_grid_ds = freq_grid_ds;
res_mfile.gmm_centroids = mu_freqspace;
res_mfile.gmm_classifications = clusts_by_chan;
res_mfile.gmm_transitions = all_transitions;

%% Try k-means on full state vector
% 
% n_comp_options = 4:8;
% AICs = zeros(1, length(n_comp_options));
% 
% for kK = 1:length(n_comp_options)
% 
%     [idx, C, sumd, D] = kmeans(U_combined, n_comp_options(kK));
%     C_freqspace = C * V_combined';
% 
%     figure;
%     sanePColor(1:n_comp_options(kK), freq_grid_ds, log(C_freqspace'), false, true);
%     set(gca, 'YScale', 'log');
%     ylabel('Frequency (Hz)');
%     xlabel('K-means centroid #');
%     title('K-means centroids (log power)');
% 
%     idx_by_chan = cell(size(all_data));
% 
%     for kC = 1:n_chans
%         idx_by_chan{kC} = nan(length(b_nan{kC}), 1);
%         idx_by_chan{kC}(~b_nan{kC}) = idx(chan_inds_nonan{kC});
%     end
% 
%     AICs(kK) = sum(sumd) + 2 * numel(C);
% end
