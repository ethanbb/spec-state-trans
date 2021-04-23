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
% score_dist_analysis(nmf_mfiles, 'L2_dist');

%% Use CCA instead
score_cca(nmf_mfiles);

%% Transition synchrony analysis

n_bootstrap = 1000;
cdf_interp_vals = 0:0.01:1;

pair_sync_scores = cell(length(exp_info), 1);
pair_sync_scores_boot = cell(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    trans_tables = cell(this_ndays, 1);
    pair_sync_scores{kE} = cell(this_ndays, 1);
    pair_sync_scores_boot{kE} = cell(this_ndays, 1);
    
    day_chan_names = cell(this_ndays, 1);
    cum_trans_interp = zeros(1, length(cdf_interp_vals));
    for kD = 1:this_ndays
        this_day = this_info.input_s(kD).name;
        mt_paths = this_info.input_s(kD).mt_res_in;
        nmf_path = this_info.input_s(kD).nmf_res_out;
        nmf_matfile = matfile(nmf_path);
        day_chan_names{kD} = nmf_matfile.chan_names;
        
        trans_tables{kD} = get_state_transitions(mt_paths, nmf_path);

        % add SPIKE-Sync score
        trans_tables{kD} = calc_transition_synchronization(trans_tables{kD});
        
        % add to unnormalized CDF
        n_chans = length(day_chan_names{kD});
        cum_edges = 0:1/n_chans:1;
        cum_vals = 0:1/(n_chans-1):1;
        cum_trans = histcounts(trans_tables{kD}.sync_score, cum_edges, 'Normalization', 'cumcount');
        cum_trans_interp = cum_trans_interp + interp1(cum_vals, cum_trans, cdf_interp_vals);

        % make individual plot
        fh = figure;
        plot_transitions(trans_tables{kD});
        title(sprintf('Transitions with sync scores - %s', this_info.input_s(kD).name));
        savefig(fh, fullfile(sr_dirs.results, this_info.days{kD}, ...
            sprintf('transitions_w_sync_%s.fig', this_info.days{kD})));
        
        % also get pairwise mean sync scores    
        pair_sync_scores{kE}{kD} = calc_pair_sync_scores(trans_tables{kD}, day_chan_names{kD});
        
        % make individual plot
        this_hr_chan_names = util.make_hr_chan_names(day_chan_names{kD}, 140);
        fh = plot_dist_mat(pair_sync_scores{kE}{kD}, this_hr_chan_names, this_day, 'trans_sync_scores');
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('pair_sync_scores_%s.fig', this_day)));
    end

    % Make CDF of sync scores and compare to bootstrap
    fh = figure;
    cdf_trans = cum_trans_interp / cum_trans_interp(end);
    plot(cdf_interp_vals, cdf_trans);
    hold on;

    cum_trans_boot = zeros(n_bootstrap, length(cdf_interp_vals));
    for kD = 1:this_ndays
        pair_sync_scores_boot{kE}{kD} = zeros([size(pair_sync_scores{kE}{kD}), n_bootstrap]);
    end
    
    for kB = 1:n_bootstrap
        shuffled_sync_scores = cell(this_ndays, 1);

        for kD = 1:this_ndays
            trans_shuffled = shuffle_trans_dwell_times(trans_tables{kD});
            trans_shuffled_w_sync = calc_transition_synchronization(trans_shuffled);
            
            % add to unnormalized CDF
            n_chans = length(day_chan_names{kD});
            cum_edges = 0:1/n_chans:1;
            cum_vals = 0:1/(n_chans-1):1;
            cum_trans = histcounts(trans_shuffled_w_sync.sync_score, cum_edges, 'Normalization', 'cumcount');
            cum_trans_boot(kB, :) = cum_trans_boot(kB, :) + interp1(cum_vals, cum_trans, cdf_interp_vals);
            
            pair_sync_scores_boot{kE}{kD}(:, :, kB) = calc_pair_sync_scores(trans_shuffled, day_chan_names{kD});
        end
    end
    cdf_trans_boot = cum_trans_boot ./ cum_trans_boot(:, end);

    % Plot median and 95% CI of bootstrap CDFs with error bars
    boot_cdf_quantiles = quantile(cdf_trans_boot, [0.025, 0.5, 0.975]);
%     neg = boot_cdf_quantiles(2,:) - boot_cdf_quantiles(1, :);
%     pos = boot_cdf_quantiles(3,:) - boot_cdf_quantiles(2,:);
%     errorbar(cdf_interp_vals, boot_cdf_quantiles(2,:), neg, pos, '.-');
    xconf = [cdf_interp_vals, cdf_interp_vals(end:-1:1)];
    yconf = [boot_cdf_quantiles(3, :), boot_cdf_quantiles(1, end:-1:1)];
    plot(cdf_interp_vals, boot_cdf_quantiles(2, :), 'r');
    fill(xconf, yconf, 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    title(sprintf('Synchronization scores for all days (%s) (%d bootstraps)', ...
        this_info.type, n_bootstrap), 'Interpreter', 'none');
    xlabel('Sync score');
    ylabel('CDF');
    legend('Real', 'Shuffled', '95% CI', 'Location', 'northwest');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
        sprintf('transition_sync_cdf_boot_%s.fig', this_info.type)));
end

%% Construct matrices of mutual information and canonical correlation, with bootsrap versions

all_chan_names = cell(length(exp_info), 1); 
hr_chan_names = cell(length(exp_info), 1); % human-readable
mut_info_combined = cell(length(exp_info), 1);
mut_info_boot = cell(length(exp_info), 1);
all_cca = cell(length(exp_info), 1);
all_cca_boot = cell(length(exp_info), 1);
% all_dists_null = cell(length(exp_info), 1);
% all_recon_err = cell(length(exp_info), 1);
all_pair_sync_scores = cell(length(exp_info), 1);
all_pair_sync_scores_boot = cell(length(exp_info), 1);
bootstrap_rstates = cell(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    nmf_res_paths = {this_info.input_s.nmf_res_out}';
    nmf_res_mfiles = cellfun(@matfile, nmf_res_paths, 'uni', false);

    all_chan_names{kE} = this_info.all_chan_names;
    hr_chan_names{kE} = util.make_hr_chan_names(all_chan_names{kE}, 140);
    n_chans = length(all_chan_names{kE});

    % gather distance information from each recording
    mut_info_combined{kE} = nan(n_chans, n_chans, this_ndays);
    mut_info_boot{kE} = nan(n_chans, n_chans, this_ndays, n_bootstrap);
    all_cca{kE} = nan(n_chans, n_chans, this_ndays);
    all_cca_boot{kE} = nan(n_chans, n_chans, this_ndays, n_bootstrap);
%     all_recon_err{kE} = nan(n_chans, n_chans, this_ndays);
    all_pair_sync_scores{kE} = nan(n_chans, n_chans, this_ndays);
    all_pair_sync_scores_boot{kE} = nan(n_chans, n_chans, this_ndays, n_bootstrap);
    bootstrap_rstates{kE} = cell(this_ndays, 1);
    
    for kD = 1:this_ndays
        this_day = this_info.input_s(kD).name;
        this_mfile = nmf_res_mfiles{kD};
        this_chans = this_mfile.all_chans;
                
        % compute mutual information
        classes = this_mfile.nmf_classes;
%         % just use run 1 of 2 (arbitrarily)
%         classes_cell = classes_cell{1};
%         classes = horzcat(classes_cell{:});
%         [~, norm_mut_info] = class_mut_info(classes);
        % actually use both runs
        % since MI is symmetric, doesn't matter which run is which argument
        classes_a = horzcat(classes{1}{:});
        classes_b = horzcat(classes{2}{:});
        [~, norm_mut_info] = class_mut_info(classes_a, classes_b);
        
        % plot mutual information
        this_hr_chan_names = util.make_hr_chan_names(this_chans, 140);
        fh = plot_dist_mat(norm_mut_info, this_hr_chan_names, this_day, 'norm_mutual_info');
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('norm_mut_info_%s.fig', this_day)));
        
        % save data to 3D array
        insert_inds = cellfun(@(c) find(strcmp(all_chan_names{kE}, c)), this_chans);
        mut_info_combined{kE}(insert_inds, insert_inds, kD) = norm_mut_info;
%         all_dists{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dists;
%         all_dists_null{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dists_null;
%         all_recon_err{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dist_recon_err;
        all_cca{kE}(insert_inds, insert_inds, kD) = this_mfile.all_cca_sim;
        all_pair_sync_scores{kE}(insert_inds, insert_inds, kD) = pair_sync_scores{kE}{kD};
        all_pair_sync_scores_boot{kE}(insert_inds, insert_inds, kD, :) = ...
            permute(pair_sync_scores_boot{kE}{kD}, [1, 2, 4, 3]);
        
        % also do bootstraps
        this_n_chans = length(this_chans);
        nmf_V = this_mfile.nmf_V;
        trans = this_mfile.nmf_transitions;
        models = cell(this_n_chans, 1);
        rstates = cell(n_bootstrap, this_n_chans);

        for kB = 1:n_bootstrap
            V_shuffled = cell(this_n_chans, 1);
            classes_shuffled = zeros(size(classes_a));

            for kC = 1:this_n_chans
                [V_shuffled{kC}, classes_shuffled(:, kC), models{kC}, rstates{kB, kC}] = ...
                    util.shuffle_scores_markov(nmf_V{1}{kC}, classes{1}{kC}, trans{1}{kC}, models{kC});
            end
            
            % normalized mutual information
            [~, mut_info_boot{kE}(insert_inds, insert_inds, kD, kB)] = ...
                class_mut_info(classes_shuffled, classes_b);
            
            % canonical correlation
            this_cca = zeros(this_n_chans);
            for iC = 1:this_n_chans
                for jC = 1:this_n_chans
                    [~, ~, ps] = canoncorr(V_shuffled{iC}, nmf_V{2}{jC});
                    this_cca(iC, jC) = mean(ps);
                end
            end
            all_cca_boot{kE}(insert_inds, insert_inds, kD, kB) = this_cca;
        end        
        bootstrap_rstates{kE}{kD} = rstates;
    end
    
    % find channels with no data
    nonempty_chan_combos = any(~isnan(mut_info_combined{kE}), 3);
    nonempty_chans = any(nonempty_chan_combos);
    
    all_chan_names{kE} = all_chan_names{kE}(nonempty_chans);
    hr_chan_names{kE} = hr_chan_names{kE}(nonempty_chans);
    mut_info_combined{kE} = mut_info_combined{kE}(nonempty_chans, nonempty_chans, :);
    mut_info_boot{kE} = mut_info_boot{kE}(nonempty_chans, nonempty_chans, :, :);
%     all_dists{kE} = all_dists{kE}(nonempty_chans, nonempty_chans, :);
%     all_dists_null{kE} = all_dists_null{kE}(nonempty_chans, nonempty_chans, :);
%     all_recon_err{kE} = all_recon_err{kE}(nonempty_chans, nonempty_chans, :);
    all_cca{kE} = all_cca{kE}(nonempty_chans, nonempty_chans, :);
    all_cca_boot{kE} = all_cca_boot{kE}(nonempty_chans, nonempty_chans, :, :);
    all_pair_sync_scores{kE} = all_pair_sync_scores{kE}(nonempty_chans, nonempty_chans, :);
    all_pair_sync_scores_boot{kE} = all_pair_sync_scores_boot{kE}(nonempty_chans, nonempty_chans, :, :);
end

%%
save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'n_bootstrap', 'all_chan_names', 'hr_chan_names', ...
    'mut_info_combined', 'mut_info_boot', 'all_cca', 'all_cca_boot', ...
    'all_pair_sync_scores', 'all_pair_sync_scores_boot', 'bootstrap_rstates', '-v7.3');

%% Mutual information combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    
    % matrix plot
    med_mut_info = median(mut_info_combined{kE}, 3, 'omitnan');
    fh = plot_dist_mat(med_mut_info, hr_chan_names{kE}, ...
        sprintf('%s, median over %d days', this_info.type, this_ndays), 'norm_mutual_info');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_nmi_%s.fig', this_info.type)));
    
    % violin plot
    fh = make_violin(mut_info_combined{kE}, all_chan_names{kE}, {
        'SameRegion', 'CrossRegion', ...
        'SameRegionL4', 'SameRegionNonL4', ...
        'CrossRegionL4', 'CrossRegionNonL4'
        }, {
        'Within region', 'Across regions', ...
        'Within, pairs with L4', 'Within, no L4', ...
        'Across, pairs with L4', 'Across, no L4'
        });
    xtickangle(20);
    ylabel('Normalized MI');
    title(sprintf('Normalized mutual information of channel classes (%s)', this_info.type), ...
        'Interpreter', 'none');
    
    pvals = ones(3, 1);
        
    % calc p-values using bootstraps which is really a permutation test, oops
    pvals(1) = permutation_tval_test(mut_info_combined{kE}, mut_info_boot{kE}, all_chan_names{kE}, ...
        'SameRegion', 'CrossRegion');
    exp_info(kE).nmi_crossvssame_pval = pvals(1);
    
    pvals(2) = permutation_tval_test(mut_info_combined{kE}, mut_info_boot{kE}, all_chan_names{kE}, ...
        'SameRegionNonL4', 'SameRegionL4');
    exp_info(kE).nmi_samel4_pval = pvals(2);
    
    pvals(3) = permutation_tval_test(mut_info_combined{kE}, mut_info_boot{kE}, all_chan_names{kE}, ...
        'CrossRegionNonL4', 'CrossRegionL4');
    exp_info(kE).nmi_crossl4_pval = pvals(3);

    % significance stars
    hs = sigstar({[1, 2], [3, 4], [5, 6]}, pvals);
    set(hs(:, 2), 'FontSize', 18);
    delete(hs(pvals >= 0.05, :));
    
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('nmi_violin_%s.fig', this_info.type)));
end

%% CCA combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    % matrix plot
    med_cca = median(all_cca{kE}, 3, 'omitnan');    
    fh = plot_dist_mat(med_cca, hr_chan_names{kE}, ...
        sprintf('%s - median over %d days', this_info.type, this_ndays), 'cca_mean_r');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_cca_r_%s.fig', this_info.type)));
    
    % violin plot
    fh = make_violin(all_cca{kE}, all_chan_names{kE}, {
        'SameRegion', 'CrossRegion', ...
        'SameRegionL4', 'SameRegionNonL4', ...
        'CrossRegionL4', 'CrossRegionNonL4'
        }, {
        'Within region', 'Across regions', ...
        'Within, pairs with L4', 'Within, no L4', ...
        'Across, pairs with L4', 'Across, no L4'
        });
    xtickangle(20);
    ylabel('Mean CCA r');
    title(sprintf('Canonical correlation of channel NMF scores (%s)', this_info.type), ...
        'Interpreter', 'none');
    
    pvals = ones(3, 1);
        
    % calc p-values using bootstraps which is really a permutation test, oops
    pvals(1) = permutation_tval_test(all_cca{kE}, all_cca_boot{kE}, all_chan_names{kE}, ...
        'SameRegion', 'CrossRegion');
    exp_info(kE).cca_crossvssame_pval = pvals(1);
    
    pvals(2) = permutation_tval_test(all_cca{kE}, all_cca_boot{kE}, all_chan_names{kE}, ...
        'SameRegionNonL4', 'SameRegionL4');
    exp_info(kE).cca_samel4_pval = pvals(2);
    
    pvals(3) = permutation_tval_test(all_cca{kE}, all_cca_boot{kE}, all_chan_names{kE}, ...
        'CrossRegionNonL4', 'CrossRegionL4');
    exp_info(kE).cca_crossl4_pval = pvals(3);

    % significance stars
    hs = sigstar({[1, 2], [3, 4], [5, 6]}, pvals);
    set(hs(:, 2), 'FontSize', 18);
    delete(hs(pvals >= 0.05, :));
    
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('cca_violin_%s.fig', this_info.type)));
end

%% Transition pair synchrony combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    % matrix plot
    med_sync = median(all_pair_sync_scores{kE}, 3, 'omitnan');    
    fh = plot_dist_mat(med_sync, hr_chan_names{kE}, ...
        sprintf('%s - median over %d days', this_info.type, this_ndays), 'trans_sync_scores');
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_trans_sync_%s.fig', this_info.type)));
    
    % violin plot
    fh = make_violin(all_pair_sync_scores{kE}, all_chan_names{kE}, {
        'SameRegion', 'CrossRegion', ...
        'SameRegionL4', 'SameRegionNonL4', ...
        'CrossRegionL4', 'CrossRegionNonL4'
        }, {
        'Within region', 'Across regions', ...
        'Within, pairs with L4', 'Within, no L4', ...
        'Across, pairs with L4', 'Across, no L4'
        });
    xtickangle(20);
    ylabel('Mean SYNC score');
    title(sprintf('Synchronization of channel class transitions (%s)', this_info.type), ...
        'Interpreter', 'none');
    
    pvals = ones(3, 1);
        
    % calc p-values using bootstraps which is really a permutation test, oops
    pvals(1) = permutation_tval_test(all_pair_sync_scores{kE}, all_pair_sync_scores_boot{kE}, ...
        all_chan_names{kE}, 'SameRegion', 'CrossRegion');
    exp_info(kE).trans_sync_crossvssame_pval = pvals(1);
    
    pvals(2) = permutation_tval_test(all_pair_sync_scores{kE}, all_pair_sync_scores_boot{kE}, ...
        all_chan_names{kE}, 'SameRegionNonL4', 'SameRegionL4');
    exp_info(kE).trans_sync_samel4_pval = pvals(2);
    
    pvals(3) = permutation_tval_test(all_pair_sync_scores{kE}, all_pair_sync_scores_boot{kE}, ...
        all_chan_names{kE}, 'CrossRegionNonL4', 'CrossRegionL4');
    exp_info(kE).trans_sync_crossl4_pval = pvals(3);

    % significance stars
    hs = sigstar({[1, 2], [3, 4], [5, 6]}, pvals);
    set(hs(:, 2), 'FontSize', 18);
    delete(hs(pvals >= 0.05, :));
    
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('trans_sync_violin_%s.fig', this_info.type)));
end

%% Finally, save
save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'exp_info', '-append');

%% Score distance combined plots
% 
% for kE = 1:length(exp_info)
%     this_info = exp_info(kE);
%     this_ndays = length(this_info.days);
% 
%     med_dists = median(all_dists{kE}, 3, 'omitnan');
%     med_dists_from_null = median(all_dists_null{kE} - all_dists{kE}, 3, 'omitnan');
%     med_recon_err = median(all_recon_err{kE}, 3, 'omitnan');
% 
%     % matrix plots
%     % Score distance
%     fh = plot_dist_mat(med_dists, hr_chan_names{kE}, ...
%         sprintf('%s - median over %d days', this_info.type, this_ndays), 'L2_dist');
%     savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_score_dist_%s.fig', this_info.type)));
%     
%     % Score dist, difference from null model
%     fh = plot_dist_mat(med_dists_from_null, hr_chan_names{kE}, ...
%         sprintf('%s (median)', this_info.type), 'L2_dist_from_null');
%     savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_score_dist_fromnull_%s.fig', this_info.type)));
%     
%     % Reconstruction error
%     fh = plot_dist_mat(med_recon_err, hr_chan_names{kE}, ...
%         sprintf('%s (median)', this_info.type), 'recon_from_L2_dist');
%     savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_l2_recon_err_%s.fig', this_info.type)));
%  
%     % Score dist, ratio of real to null model
%     dist_ratio = all_dists{kE} ./ all_dists_null{kE};
%     [fh, pval] = make_crosss_vs_same_violin(dist_ratio, all_chan_names{kE}, 'right');
%     ylabel('Score distance (frac. of null model)');
%     title({'Euclidean distance of aligned NMF scores, compared to null model', ...
%         this_info.type}, 'Interpreter', 'none');
%     savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
%         sprintf('dist_vs_null_violin_%s.fig', this_info.type)));
%     exp_info(kE).dist_vs_null_pval = pval;
%     fprintf('Score distance vs. null, %s, cross- vs. same-region:\n p = %g\n', this_info.type, pval);
%     
%     % Reconstruction error
%     [fh, pval] = make_cross_vs_same_violin(all_recon_err{kE}, all_chan_names{kE}, 'right');
%     ylabel('Reconstruction error');
%     title({'Cross-channel reconstruction error from aligned NMF scores', ...
%         this_info.type}, 'Interpreter', 'none');
%     savefig(fh, fullfile(sr_dirs.results, 'res_figs', ...
%         sprintf('recon_err_violin_%s.fig', this_info.type)));
%     exp_info(kE).recon_err_pval = pval;
%     fprintf('Reconstruction error, %s, cross- vs. same-region:\n p = %g\n', this_info.type, pval);
% end

%% Function to get pairwise mean sync scores
function pair_sync_scores = calc_pair_sync_scores(trans_table, chan_names)
n_chans = length(chan_names);
pair_sync_scores = zeros(n_chans);

% get indices of each channel's transitions
chan_trans_inds = cellfun(@(cn) strcmp(trans_table.chan_name, cn), chan_names, 'uni', false);
for iC = 1:n_chans
    for jC = 1:iC-1
        % compute sync score for just transitions in one of the two current channels
        pair_trans_table = trans_table(chan_trans_inds{iC} | chan_trans_inds{jC}, :);
        pair_trans_table = calc_transition_synchronization(pair_trans_table);
        pair_sync_scores(iC, jC) = mean(pair_trans_table.sync_score);
    end
end

% the matrix should be symmetric & diagonal entries are meaningless
pair_sync_scores = pair_sync_scores + pair_sync_scores';
pair_sync_scores(linspace(1, n_chans^2, n_chans)) = nan;

end

%% Functions to make violin plots

function vals_bytype = categorize_pair_data(mats_over_days, chan_names)
% Take some matrices of data over days (in 3rd dimension) and return struct w/ combined
% data vectors, with nans removed, in these mutually exclusive categories:
% 'SameChannel', 'SameRegionL4', 'SameRegionNonL4', 'CrossRegionL4', 'CrossRegionNonL4'
% Any pair of a channel with itself is 'SameChannel'
% Any pair where either channel is labeled as L4 is in the corresponding L4 category
% Other pairs are in the corresponding NonL4 category depending on whether the region is the same.

chan_names = string(chan_names(:));
n_chans = length(chan_names);
chan_name_parts = split(chan_names, "_");
regions = chan_name_parts(:, 1);
layers = chan_name_parts(:, 2);

layer_isl4 = layers == "L4";
pair_hasl4 = layer_isl4 | layer_isl4';
same_region = regions == regions';
same_channel = logical(eye(n_chans));

% to broadcast along later dimensions
broadcast_dims = size(mats_over_days, 3:ndims(mats_over_days));
get_masked_vals = @(mask) mats_over_days(repmat(mask, [1, 1, broadcast_dims]));

vals_bytype.SameChannel = get_masked_vals(same_channel);
vals_bytype.SameRegionL4 = get_masked_vals(~same_channel & pair_hasl4 & same_region);
vals_bytype.SameRegionNonL4 = get_masked_vals(~same_channel & ~pair_hasl4 & same_region);
vals_bytype.CrossRegionL4 = get_masked_vals(~same_channel & pair_hasl4 & ~same_region);
vals_bytype.CrossRegionNonL4 = get_masked_vals(~same_channel & ~pair_hasl4 & ~same_region);
vals_bytype.SameRegion = [vals_bytype.SameRegionL4; vals_bytype.SameRegionNonL4];
vals_bytype.CrossRegion = [vals_bytype.CrossRegionL4; vals_bytype.CrossRegionNonL4];

vals_bytype = structfun(@(v) v(~isnan(v)), vals_bytype, 'uni', false);

end

function [hf, vs] = make_violin(mats_over_days, chan_names, category_snames, category_hr_names)
% sname = struct name, hr_name = human-readable name

violin_colors = struct(...
    'SameChannel', [0.5 0.5 0.5], ...
    'SameRegion', [0 0 1], ...
    'SameRegionL4', [0.5 0.5 1], ...
    'SameRegionNonL4', [0 0 0.5], ...
    'CrossRegion', [1 0 0], ...
    'CrossRegionL4', [1 0.5 0.5], ...
    'CrossRegionNonL4', [0.5 0 0] ...
    );

assert(length(category_snames) == length(category_hr_names), 'Mismatch in # of categories');
nplots = length(category_snames);

vals_bytype = categorize_pair_data(mats_over_days, chan_names);
violin_s = struct;
my_colors = zeros(nplots, 3);
for kV = 1:nplots
    sname = category_snames{kV};
    violin_s.(sname) = vals_bytype.(sname);
    my_colors(kV, :) = violin_colors.(sname);
end

hf = figure;
vs = violinplot(violin_s);
xticklabels(category_hr_names);

for kV = 1:nplots
    vs(kV).ViolinColor = my_colors(kV, :);
end

end

function [hf, p] = make_cross_vs_same_violin(mats_over_days, chan_names, ...
    ttest_tail, b_include_self)
% make violin plot to compare a quantity for cross-region pairs, witin-region pairs, and self-loops.
% chan_names should be in "region_layer" format, not the "human-readable" format.
% ttest_tail: 'right' means hypothesis is across > within for this measure and vice versa.

if ~exist('ttest_tail', 'var') || isempty(ttest_tail)
    ttest_tail = 'left';
end

if ~exist('b_include_self', 'var') || isempty(b_include_self)
    b_include_self = true;
end

types = {'Single channel replication', 'Within region', 'Across regions'};
type_snames = {'SameChannel', 'SameRegion', 'CrossRegion'};
if ~b_include_self
    types = types(2:end);
    type_snames = type_snames(2:end);
end

[hf, vs] = make_violin(mats_over_days, chan_names, type_snames, types);

if b_include_self
    same_region_data = vs(2).ScatterPlot.YData;
    cross_region_data = vs(3).ScatterPlot.YData;
else
    same_region_data = vs(1).ScatterPlot.YData;
    cross_region_data = vs(2).ScatterPlot.YData;
end
    
[~, p] = ttest2(cross_region_data, same_region_data, 0.01, ttest_tail);

end

function p = permutation_tval_test(real_mats, perm_mats, chan_names, contrast_group1, contrast_group2)
% Get permutation test p-value for the hypothesis contrast_group1 > contrast_group2
% These are strings that refer to struct names of the return value of categorize_pair_data
% real_mats should be a 3-dimensional array of the measure on channel pairs across days
% perm_mats should be a 4-dimensional array, the last dimension iterates over permutation runs
% chan_names is the machine-readable channel names to input to categorize_pair_data

welch_tval = @(group1, group2) (mean(group1) - mean(group2)) / ...
    sqrt(var(group1) / length(group1) + var(group2) / length(group2));
tval_fn = @(vals_bytype) welch_tval(vals_bytype.(contrast_group1), vals_bytype.(contrast_group2));

% get real tval
real_vals_bytype = categorize_pair_data(real_mats, chan_names);
tval_real = tval_fn(real_vals_bytype);

% get permutation tvals
perm_vals_bytype = cellfun(@(perm_cca) categorize_pair_data(perm_cca, chan_names), num2cell(perm_mats, [1, 2, 3]));
tvals_perm = arrayfun(tval_fn, perm_vals_bytype);

p = sum(tvals_perm > tval_real) / size(perm_mats, 4);

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

