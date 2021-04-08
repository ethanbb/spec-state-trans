% Cluster states for each channel on each day using non-negative matrix factorization, then
% generate null model data, find KL divergence between pairs of channels (after transformation),
% and make plots.

sr_dirs = prepSR;

% set up full channel name lists for analyses later
layer_names = [
    arrayfun(@(k) strcat("Sup", num2str(k)), (8:-1:1)')
    "L4"
    arrayfun(@(k) strcat("Inf", num2str(k)), (1:8)')
];

make_layer_names = @(regions) reshape(string(regions(:)') + "_" + layer_names, [], 1);

exp_info = [
    struct('type', 'm1_v1', 'all_chan_names', make_layer_names(["M1", "V1"]), ...
        'days', {{
%             '2020-01-30'
%             '2020-01-31'
%             '2020-02-06'
%             '2020-03-05'
%             '2020-03-06'
%             '2020-03-10'
%             '2020-03-11'
            % not using any of the above - short probe recordings
            '2020-10-26'
            '2020-10-27'
            '2020-10-28'
            '2020-10-29'
        }})
    struct('type', 'bilateral', 'all_chan_names', make_layer_names(["V1L", "V1R"]), ...
        'days', {{
            '2021-01-27'
            '2021-01-29'
            '2021-01-31'
            '2021-02-02'
        }})
    struct('type', 'm1_v1_csd', 'all_chan_names', make_layer_names(["M1", "V1"]), ...
        'days', {{
            '2020-10-26_csd'
            '2020-10-27_csd'
            '2020-10-28_csd'
            '2020-10-29_csd'
        }})
    struct('type', 'bilateral_csd', 'all_chan_names', make_layer_names(["V1L", "V1R"]), ...
        'days', {{
            '2021-01-27_csd'
            '2021-01-29_csd'
            '2021-01-31_csd'
            '2021-02-02_csd'
        }})
    ];

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

%% Euclidean distance analysis of aligned scores
score_dist_analysis(nmf_mfiles, 'L2_dist');

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

%% Construct matrices of mutual information and score distance

chan_names = cell(length(exp_info), 1); 
hr_chan_names = cell(length(exp_info), 1); % human-readable
mut_info_combined = cell(length(exp_info), 1);
all_dists = cell(length(exp_info), 1);
all_dists_null = cell(length(exp_info), 1);
all_recon_err = cell(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    nmf_res_paths = {this_info.input_s.nmf_res_out}';
    nmf_res_mfiles = cellfun(@matfile, nmf_res_paths, 'uni', false);

    chan_names{kE} = this_info.all_chan_names;
    hr_chan_names{kE} = util.make_hr_chan_names(chan_names{kE}, 140);
    n_chans = length(chan_names{kE});

    % gather distance information from each recording
    mut_info_combined{kE} = nan(n_chans, n_chans, this_ndays);
    all_dists{kE} = nan(n_chans, n_chans, this_ndays);
    all_dists_null{kE} = nan(n_chans, n_chans, this_ndays);
    all_recon_err{kE} = nan(n_chans, n_chans, this_ndays);
    
    for kD = 1:this_ndays
        this_day = this_info.input_s(kD).name;
        this_mfile = nmf_res_mfiles{kD};
        this_chans = this_mfile.all_chans;
                
        % compute mutual information
        classes_cell = this_mfile.nmf_classes;
        % just use run 1 of 2 (arbitrarily)
        classes_cell = classes_cell{1};
        classes = horzcat(classes_cell{:});
        [~, norm_mut_info] = class_mut_info(classes);
        
        % plot mutual information
        this_hr_chan_names = util.make_hr_chan_names(this_chans, 140);
        fh = plot_dist_mat(norm_mut_info, this_hr_chan_names, this_day, 'norm_mutual_info');
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('norm_mut_info_%s.fig', this_day)));
        
        % save data to 3D array
        insert_inds = cellfun(@(c) find(strcmp(chan_names{kE}, c)), this_chans);
        mut_info_combined{kE}(insert_inds, insert_inds, kD) = norm_mut_info;
        all_dists{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dists;
        all_dists_null{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dists_null;
        all_recon_err{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dist_recon_err;
    end
    
    % find channels with no data
    nonempty_chan_combos = any(~isnan(mut_info_combined{kE}), 3);
    nonempty_chans = any(nonempty_chan_combos);
    
    chan_names{kE} = chan_names{kE}(nonempty_chans);
    hr_chan_names{kE} = hr_chan_names{kE}(nonempty_chans);
    mut_info_combined{kE} = mut_info_combined{kE}(nonempty_chans, nonempty_chans, :);
    all_dists{kE} = all_dists{kE}(nonempty_chans, nonempty_chans, :);
    all_dists_null{kE} = all_dists_null{kE}(nonempty_chans, nonempty_chans, :);
    all_recon_err{kE} = all_recon_err{kE}(nonempty_chans, nonempty_chans, :);
end

%% Mutual information combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    
    med_mut_info = median(mut_info_combined{kE}, 3, 'omitnan');

    % matrix plot
    fh = plot_dist_mat(med_mut_info, hr_chan_names{kE}, ...
        sprintf('%s, median of %d days', this_info.type, this_ndays), 'norm_mutual_info');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('norm_mut_info_%s.fig', this_info.type)));
    
    % violin plot
    [fh, pval] = make_violin_from_mats(mut_info_combined{kE}, chan_names{kE}, 'left', false);
    ylabel('Normalized MI');
    title({'Normalized mutual information of channel pairs', this_info.type}, 'Interpreter', 'none');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
        sprintf('norm_mut_info_violin_%s.fig', this_info.type)));
    
    % stats    
    exp_info(kE).nmf_pval = pval;
    fprintf('NMI, %s, cross- vs. same-region:\n p = %g\n', this_info.type, pval);
end

%% Score distance combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    med_dists = median(all_dists{kE}, 3, 'omitnan');
    med_dists_from_null = median(all_dists_null{kE} - all_dists{kE}, 3, 'omitnan');
    med_recon_err = median(all_recon_err{kE}, 3, 'omitnan');

    % matrix plots
    % Score distance
    fh = plot_dist_mat(med_dists, hr_chan_names{kE}, ...
        sprintf('%s - median over %d days', this_info.type, this_ndays), 'L2_dist');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_score_dist_%s.fig', this_info.type)));
    
    % Score dist, difference from null model
    fh = plot_dist_mat(med_dists_from_null, hr_chan_names{kE}, ...
        sprintf('%s (median)', this_info.type), 'L2_dist_from_null');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_score_dist_fromnull_%s.fig', this_info.type)));
    
    % Reconstruction error
    fh = plot_dist_mat(med_recon_err, hr_chan_names{kE}, ...
        sprintf('%s (median)', this_info.type), 'recon_from_L2_dist');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_l2_recon_err_%s.fig', this_info.type)));
    
    % violin plots
    % Score dist, ratio of real to null model
    dist_ratio = all_dists{kE} ./ all_dists_null{kE};
    [fh, pval] = make_violin_from_mats(dist_ratio, chan_names{kE}, 'right');
    ylabel('Score distance (frac. of null model)');
    title({'Euclidean distance of aligned NMF scores, compared to null model', ...
        this_info.type}, 'Interpreter', 'none');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
        sprintf('dist_vs_null_violin_%s.fig', this_info.type)));
    exp_info(kE).dist_vs_null_pval = pval;
    fprintf('Score distance vs. null, %s, cross- vs. same-region:\n p = %g\n', this_info.type, pval);
    
    % Reconstruction error
    [fh, pval] = make_violin_from_mats(all_recon_err{kE}, chan_names{kE}, 'right');
    ylabel('Reconstruction error');
    title({'Cross-channel reconstruction error from aligned NMF scores', ...
        this_info.type}, 'Interpreter', 'none');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
        sprintf('recon_err_violin_%s.fig', this_info.type)));
    exp_info(kE).recon_err_pval = pval;
    fprintf('Reconstruction error, %s, cross- vs. same-region:\n p = %g\n', this_info.type, pval);
end

%% Finally, save
save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'exp_info');

%% Function to make violin plots

function [hf, p] = make_violin_from_mats(mats_over_days, chan_names, ttest_tail, b_include_self)
% make violin plot to compare a quantity for cross-region pairs, witin-region pairs, and self-loops.
% chan_names should be in "region_layer" format, not the "human-readable" format.
% ttest_tail: 'right' means hypothesis is across > within for this measure and vice versa.

if ~exist('b_include_self', 'var')
    b_include_self = true;
end

types = {'Single channel replication', 'Within region', 'Across regions'};

if ~b_include_self
    types = types(2:end);
end

type_snames = cellfun(@matlab.lang.makeValidName, types, 'uni', false);
n_types = length(types);

n_chans = length(chan_names);
this_regions = strtok(chan_names, '_');
same_region = this_regions == this_regions';

% get indices of each type of pair
self_loops = find(eye(n_chans));
linear_inds = {
	setdiff(find(same_region), self_loops)
    find(~same_region)
    };

if b_include_self
    linear_inds = [{self_loops}; linear_inds];
end

vals_bytype = struct;

for kT = 1:n_types
    tname = type_snames{kT};
    
    vals_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(mats_over_days, [1, 2]), 'uni', false);
    vals_bytype.(tname) = vertcat(vals_perday{:});
end

hf = figure;
violinplot(vals_bytype);
xticklabels(types);

[~, p] = ttest2(vals_bytype.AcrossRegions, vals_bytype.WithinRegion, 0.01, ttest_tail);
end

%% old code for separate real and null dist violin plots
%     hf = figure; hold on;
%     viol1 = violinplot(dists_bytype, [], 'ViolinColor', [1, 0, 0]);
%     viol2 = violinplot(dists_null_bytype, [], 'ViolinColor', [0.5, 0.5, 0.5]);
%     legend([viol1(1).ViolinPlot, viol2(1).ViolinPlot], 'Real data', 'Markov-chain generated', 'Location', 'southeast');
%     xticklabels(types);
%     ylabel('Euclidean score distance');
%
%     title({'Euclidean distance of aligned NMF scores between pairs of channels', ...
%         exp_info(kE).type}, 'Interpreter', 'none');
%
%     savefig(hf, fullfile(sr_dirs.results, 'res_figs', ...
%         sprintf('score_dist_violin_%s.fig', exp_info(kE).type)));

%% make a graph plot out of it
% med_kl_div_from_null_nonnan = med_dists_from_null;
% med_kl_div_from_null_nonnan(isnan(med_dists_from_null)) = 0;
% g = digraph(med_dists_from_null, 'omitselfloops');
% ecolors = 5 * g.Edges.Weight / max(g.Edges.Weight);
% edge_labels = arrayfun(@(w) sprintf('%.2f', w), g.Edges.Weight, 'uni', false);
%
% hf_graph = figure;
% plot(g, 'EdgeCData', ecolors, 'EdgeColor', 'flat', 'NodeLabel', chans, ... 'EdgeLabel', edge_labels, ...
%      'Interpreter', 'none', 'NodeFontSize', 12, 'NodeFontWeight', 'bold', 'Layout', 'layered', ...
%      'Direction', 'right', 'Sources', 1:length(chans)/2, 'Sinks', length(chans)/2+1:length(chans));
%
% title(sprintf('Synchrony of channel pairs\n(bits below null model KL divergence)'));
% colorbar;
% savefig(hf_graph, fullfile(sr_dirs.results, 'res_figs', 'med_kl_div_fromnull_graphplot.fig'));

