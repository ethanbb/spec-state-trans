% Cluster states for each channel on each day using non-negative matrix factorization, then
% generate null model data, find KL divergence between pairs of channels (after transformation),
% and make plots.

sr_dirs = prepSR;

exp_info = [struct('type', 'm1_v1', 'days', {{
%         '2020-01-30'
%         '2020-01-31'
%         '2020-02-06'
%         '2020-03-05'
%         '2020-03-06'
%         '2020-03-10'
%         '2020-03-11'
        % not using any of the above - short probe recordings
        '2020-10-26'
        '2020-10-27'
        '2020-10-28'
        '2020-10-29'
    }})
    struct('type', 'bilateral', 'days', {{
        '2021-01-27'
        '2021-01-29'
        '2021-01-31'
        '2021-02-02'
    }})];

% set up full channel name lists for analyses later
layer_names = [
    arrayfun(@(k) ['Sup', num2str(k)], 8:-1:1, 'uni', false), {'L4'}, ...
    arrayfun(@(k) ['Inf', num2str(k)], 1:8, 'uni', false)
    ];

exp_info(1).all_chan_names = [strcat('M1_', layer_names), strcat('V1_', layer_names)];
exp_info(2).all_chan_names = [strcat('V1L_', layer_names), strcat('V1R_', layer_names)];

exp_types = {exp_info.type};
all_days = vertcat(exp_info.days);

n_days = length(all_days);

%% Collect datasets on each day

for kE = 1:length(exp_types)
    this_days = exp_info(kE).days;
    this_ndays = length(this_days);
    exp_info(kE).input_s =  struct(...
        'name', this_days, ...
        'mt_res_in', cell(this_ndays, 1), ...
        'nmf_res_out', cell(this_ndays, 1), ...
        'xval_fig_dir', cell(this_ndays, 1));

    for kD = 1:this_ndays
        % Build input to concat_and_nmf

        curr_day = this_days{kD};
        exp_info(kE).input_s(kD).nmf_res_out = fullfile(sr_dirs.results, curr_day, 'nmf_res.mat');
        exp_info(kE).input_s(kD).xval_fig_dir = fullfile(sr_dirs.results, curr_day);

        % find all "layers" results files under this day
        res_fn = 'mt_res_layers.mat';
        res_entries = dir(fullfile(sr_dirs.results, curr_day, '*', res_fn));
        dirs = sort({res_entries.folder});
        exp_info(kE).input_s(kD).mt_res_in = fullfile(dirs, res_fn);
    end
end

input_s_all = vertcat(exp_info.input_s);

%% Do NMF
nmf_mfiles = concat_and_nmf(input_s_all);

%% Make null model data
gen_null_model_data(nmf_mfiles);

%% Do KL divergence analysis
% score_dist_analysis(nmf_mfiles, 'kl_div');

%% Euclidean distance analysis of aligned scores
score_dist_analysis(nmf_mfiles, 'L2_dist');

%% Mutual information analysis - now normalized

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    chan_names = this_info.all_chan_names;
    n_chans = length(chan_names);
    mut_info_combined = nan(n_chans, n_chans, this_ndays);
    
    for kD = 1:this_ndays
        this_day = this_info.input_s(kD).name;
        mfile = matfile(this_info.input_s(kD).nmf_res_out);
        
        classes_cell = mfile.nmf_classes;
        % just use run 1 of 2 (arbitrarily)
        classes_cell = classes_cell{1};
        chans = mfile.all_chans;

        classes = horzcat(classes_cell{:});
        [~, norm_mut_info] = class_mut_info(classes);

        fh = plot_dist_mat(norm_mut_info, chans, this_day, 'norm_mutual_info');
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('norm_mut_info_%s.fig', this_day)));
        
        % save to combined array
        insert_inds = cellfun(@(c) find(strcmp(chan_names, c)), chans);
        mut_info_combined(insert_inds, insert_inds, kD) = norm_mut_info;
    end
    
    % Take median and make combined plot
    med_mut_info = nanmedian(mut_info_combined, 3);
    
    % eliminate channels with no data
    chans_empty = all(isnan(med_mut_info)) & all(isnan(med_mut_info'));
    chan_names = chan_names(~chans_empty);
    med_mut_info = med_mut_info(~chans_empty, ~chans_empty);
    
    fh = plot_dist_mat(med_mut_info, chan_names, ...
        sprintf('%s, median of %d days', this_info.type, this_ndays), 'norm_mutual_info');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('norm_mut_info_%s.fig', this_info.type)));
end

%% Transition synchrony analysis
for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    trans_tables = cell(this_ndays, 1);
    
    for kD = 1:this_ndays
        mt_paths = this_info.input_s(kD).mt_res_in;
        nmf_path = this_info.input_s(kD).nmf_res_out;
        trans_tables{kD} = get_state_transitions(mt_paths, nmf_path);

        % add SPIKE-Sync score
        trans_tables{kD} = calc_transition_synchronization(trans_tables{kD});

        fh = figure;
        plot_transitions(trans_tables{kD});
        title(sprintf('Transitions with sync scores - %s', this_info.input_s(kD).name));
        savefig(fh, fullfile(sr_dirs.results, this_info.days{kD}, ...
            sprintf('transitions_w_sync_%s.fig', this_info.days{kD})));
    end

    % Make CDF of sync scores and compare to bootstrap
    n_bootstrap = 1000;

    fh = figure;
    sync_scores_all = cellfun(@(tt) tt.sync_score, trans_tables, 'uni', false);
    sync_scores_all = vertcat(sync_scores_all{:});
    [real_cdf, bin_edges] = histcounts(sync_scores_all, 'Normalization', 'cdf');
    bin_centers = mean([bin_edges(1:end-1); bin_edges(2:end)]);
    plot(bin_centers, real_cdf, '.-');
    hold on;

    boot_cdf = zeros(n_bootstrap, length(bin_centers));
    for kB = 1:n_bootstrap
        shuffled_sync_scores = cell(this_ndays, 1);

        for kD = 1:this_ndays
            trans_shuffled = shuffle_trans_dwell_times(trans_tables{kD});
            trans_shuffled = calc_transition_synchronization(trans_shuffled);
            shuffled_sync_scores{kD} = trans_shuffled.sync_score;
        end
        shuffled_sync_scores = vertcat(shuffled_sync_scores{:});
        boot_cdf(kB, :) = histcounts(shuffled_sync_scores, bin_edges, 'Normalization', 'cdf');
    end

    % Plot median and 95% CI of bootstrap CDFs with error bars
    boot_cdf_quantiles = quantile(boot_cdf, [0.025, 0.5, 0.975]);
    neg = boot_cdf_quantiles(2,:) - boot_cdf_quantiles(1, :);
    pos = boot_cdf_quantiles(3,:) - boot_cdf_quantiles(2,:);
    errorbar(bin_centers, boot_cdf_quantiles(2,:), neg, pos, '.-');
    title(sprintf('Synchronization scores for all days (%s) (%d bootstraps)', ...
        this_info.type, n_bootstrap), 'Interpreter', 'none');
    xlabel('Sync score');
    ylabel('CDF');
    legend('Real', 'Shuffled', 'Location', 'northwest');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
        sprintf('transition_sync_cdf_boot_%s.fig', this_info.type)));
end

%% Get and plot median score distance for each channel pair

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    nmf_res_paths = {this_info.input_s.nmf_res_out}';
    nmf_res_mfiles = cellfun(@matfile, nmf_res_paths, 'uni', false);

    chan_names = this_info.all_chan_names;
    n_chans = length(chan_names);

    % gather distance information from each recording
    all_dists = nan(n_chans, n_chans, this_ndays);
    all_dists_null = nan(n_chans, n_chans, this_ndays);

    for kD = 1:this_ndays
        insert_inds = cellfun(@(c) find(strcmp(chan_names, c)), nmf_res_mfiles{kD}.chan_names);

        this_mfile = nmf_res_mfiles{kD};
        all_dists(insert_inds, insert_inds, kD) = this_mfile.allL2_dists;
        all_dists_null(insert_inds, insert_inds, kD) = this_mfile.allL2_dists_null;
    end

    med_dists = nanmedian(all_dists, 3);

    % eliminate channels with no data
    chans_empty = all(isnan(med_dists)) & all(isnan(med_dists'));
    chan_names = chan_names(~chans_empty);
    med_dists = med_dists(~chans_empty, ~chans_empty);
    all_dists = all_dists(~chans_empty, ~chans_empty, :);
    all_dists_null = all_dists_null(~chans_empty, ~chans_empty, :);

    hf = plot_dist_mat(med_dists, chan_names, ...
        sprintf('%s - median over %d days', replace(this_info.type, '_', '\_'), this_ndays), 'L2_dist');
    savefig(hf, fullfile(sr_dirs.results, 'res_figs', sprintf('med_kl_div_%s.fig', this_info.type)));
    
    % Same thing but with difference from null model
    med_dists_from_null = nanmedian(all_dists_null - all_dists, 3);
    hf = plot_dist_mat(med_dists_from_null, chan_names, ...
        sprintf('%s (median)', replace(this_info.type, '_', '\_')), 'L2_dist_from_null');
    savefig(hf, fullfile(sr_dirs.results, 'res_figs', sprintf('med_kl_div_fromnull_%s.fig', this_info.type)));
end

%% **** Below here is old - to be updated ****

%% make a graph plot out of it
med_kl_div_from_null_nonnan = med_dists_from_null;
med_kl_div_from_null_nonnan(isnan(med_dists_from_null)) = 0;
g = digraph(med_dists_from_null, 'omitselfloops');
ecolors = 5 * g.Edges.Weight / max(g.Edges.Weight);
edge_labels = arrayfun(@(w) sprintf('%.2f', w), g.Edges.Weight, 'uni', false);

hf_graph = figure;
plot(g, 'EdgeCData', ecolors, 'EdgeColor', 'flat', 'NodeLabel', chans, ... 'EdgeLabel', edge_labels, ...
     'Interpreter', 'none', 'NodeFontSize', 12, 'NodeFontWeight', 'bold', 'Layout', 'layered', ...
     'Direction', 'right', 'Sources', 1:length(chans)/2, 'Sinks', length(chans)/2+1:length(chans));

title(sprintf('Synchrony of channel pairs\n(bits below null model KL divergence)'));
colorbar;
savefig(hf_graph, fullfile(sr_dirs.results, 'res_figs', 'med_kl_div_fromnull_graphplot.fig'));

%% Make violin plot

types = {'Same channel', 'Same region', 'Cross-region'};
type_snames = cellfun(@matlab.lang.makeValidName, types, 'uni', false);
n_types = length(types);

linear_inds = cell(n_types, 1);
linear_inds{1} = find(eye(n_chans));
linear_inds{2} = find(blkdiag(ones(n_chans/2), ones(n_chans/2)) - eye(n_chans));
linear_inds{end} = setdiff(1:n_chans^2, vertcat(linear_inds{1:end-1}));

kl_divs_bytype = struct;
kl_divs_null_bytype = struct;

for kT = 1:n_types
    %mean_kldiv_oftype = @(mat) mean(mat(linear_inds{kT}));
    
    kl_divs_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_dists, [1, 2]), 'uni', false);
    kl_divs_bytype.(type_snames{kT}) = vertcat(kl_divs_perday{:});
    
    kl_divs_null_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_dists_null, [1, 2]), 'uni', false);
    kl_divs_null_bytype.(type_snames{kT}) = vertcat(kl_divs_null_perday{:});
end

hf = figure; hold on;
viol1 = violinplot(kl_divs_bytype, [], 'ViolinColor', [1, 0, 0]);
viol2 = violinplot(kl_divs_null_bytype, [], 'ViolinColor', [0.5, 0.5, 0.5]);
legend([viol1(1).ViolinPlot, viol2(1).ViolinPlot], 'Real data', 'Markov-chain generated', 'Location', 'southeast');
xticklabels(types);
ylabel('KL divergence (bits)');

title({'Mean KL divergence of aligned NMF scores between pairs of channels,', ...
    'with second channel either real data or resampled based on discrete Markov chain'});

savefig(hf, fullfile(sr_dirs.results, 'res_figs', 'kl_div_withnull_violin.fig'));
