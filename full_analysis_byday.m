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

%% Make null model data (currently doing bootstrap runs later instead)
% gen_null_model_data(nmf_mfiles);

%% Euclidean distance analysis of aligned scores
% score_dist_analysis(nmf_mfiles, 'L2_dist');

%% Use CCA instead
score_cca(nmf_mfiles);

%% Construct matrices of mutual information and canonical correlation, with bootsrap versions

n_bootstrap = 1000;
cdf_interp_vals = 0:0.01:1;

all_chan_names = cell(length(exp_info), 1); 
hr_chan_names = cell(length(exp_info), 1); % human-readable
mut_info_combined = cell(length(exp_info), 1);
mut_info_boot = cell(length(exp_info), 1);
all_cca = cell(length(exp_info), 1);
all_cca_boot = cell(length(exp_info), 1);
all_pair_sync_scores = cell(length(exp_info), 1);
all_pair_sync_scores_boot = cell(length(exp_info), 1);
bootstrap_rstates = cell(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    all_chan_names{kE} = this_info.all_chan_names;
    hr_chan_names{kE} = util.make_hr_chan_names(all_chan_names{kE}, 140);
    n_chans = length(all_chan_names{kE});

    % gather distance information from each recording
    mut_info_combined{kE} = nan(n_chans, n_chans, this_ndays);
    mut_info_boot{kE} = nan(n_chans, n_chans, this_ndays, n_bootstrap);
    all_cca{kE} = nan(n_chans, n_chans, this_ndays);
    all_cca_boot{kE} = nan(n_chans, n_chans, this_ndays, n_bootstrap);
    all_pair_sync_scores{kE} = nan(n_chans, n_chans, this_ndays);
    all_pair_sync_scores_boot{kE} = nan(n_chans, n_chans, this_ndays, n_bootstrap);
    bootstrap_rstates{kE} = cell(this_ndays, 1);
    
    cum_trans_interp = zeros(1, length(cdf_interp_vals));
    cum_trans_boot = zeros(n_bootstrap, length(cdf_interp_vals));
    
    for kD = 1:this_ndays
        this_day = this_info.input_s(kD).name;
        mt_paths = this_info.input_s(kD).mt_res_in;
        this_mfile = matfile(this_info.input_s(kD).nmf_res_out, 'Writable', true);
        this_chans = this_mfile.all_chans;
        this_n_chans = length(this_chans);
                
        % compute transition sync scores
        trans_table = get_state_transitions(mt_paths, this_mfile);
        trans_table = calc_transition_synchronization(trans_table);
        
        % add to unnormalized CDF
        cum_edges = 0:1/this_n_chans:1;
        cum_vals = 0:1/(this_n_chans-1):1;
        cum_trans = histcounts(trans_table.sync_score, cum_edges, 'Normalization', 'cumcount');
        cum_trans_interp = cum_trans_interp + interp1(cum_vals, cum_trans, cdf_interp_vals);
        
        % make individual transition plot
        fh = figure;
        plot_transitions(trans_table);
        title(sprintf('Transitions with sync scores - %s', this_day), 'Interpreter', 'none');
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('transitions_w_sync_%s.fig', this_day)));
        
        % also get pairwise mean scores
        pair_sync_scores = calc_pair_sync_scores(trans_table, this_chans);
        
        % make individual pairwise sync score plot
        this_hr_chan_names = util.make_hr_chan_names(this_chans, 140);
        fh = plot_dist_mat(pair_sync_scores, this_hr_chan_names, this_day, ...
            'trans_sync_scores', 'lower_nodiag', this_chans);
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('pair_sync_scores_%s.fig', this_day)));
        
        % compute mutual information
        classes = this_mfile.filtered_classes;

        % use 1st run
        classes = horzcat(classes{1}{:});
        [~, norm_mut_info] = class_mut_info(classes);
        % since this is symmetric and we don't want to duplicate samples,
        % nan-out the upper triangular part
        norm_mut_info(1:this_n_chans >= (1:this_n_chans)') = nan;
        
        % plot mutual information
        fh = plot_dist_mat(norm_mut_info, this_hr_chan_names, this_day, 'norm_mutual_info', ...
            'lower_nodiag', this_chans);
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('norm_mut_info_%s.fig', this_day)));
        
        % save data to 3D array
        insert_inds = cellfun(@(c) find(strcmp(all_chan_names{kE}, c)), this_chans);
        mut_info_combined{kE}(insert_inds, insert_inds, kD) = norm_mut_info;
        all_cca{kE}(insert_inds, insert_inds, kD) = this_mfile.all_cca_sim;
        all_pair_sync_scores{kE}(insert_inds, insert_inds, kD) = pair_sync_scores;
        
        % also do bootstraps
        nmf_V = this_mfile.nmf_V;
        trans = this_mfile.filtered_transitions;
        models = cell(this_n_chans, 1);
        rstates = cell(n_bootstrap, this_n_chans);

        for kB = 1:n_bootstrap
            V_shuffled = cell(this_n_chans, 1);
            classes_shuffled = cell(this_n_chans, 1);

            for kC = 1:this_n_chans
                [V_shuffled{kC}, classes_shuffled{kC}, models{kC}, rstates{kB, kC}] = ...
                    util.shuffle_scores_markov(nmf_V{1}{kC}, classes(:, kC), trans{1}{kC}, models{kC});
            end
            
            % mean transition synchrony CDF
            % make struct to stand in for NMF mfile in get_state_transitions
            boot_nmf_info = struct;
            boot_nmf_info.nmf_classes = {classes_shuffled};
            boot_nmf_info.nmf_V = {V_shuffled};
            boot_nmf_info.time_axis = this_mfile.time_axis;
            boot_nmf_info.chan_names = this_mfile.chan_names;
            boot_trans_table = get_state_transitions(mt_paths, boot_nmf_info, ...
                struct('save_filtered_classes', false));
            boot_trans_table = calc_transition_synchronization(boot_trans_table);
            
            cum_trans = histcounts(boot_trans_table.sync_score, cum_edges, 'Normalization', 'cumcount');
            cum_trans_boot(kB, :) = cum_trans_boot(kB, :) + interp1(cum_vals, cum_trans, cdf_interp_vals);
            
            % pair sync scores
            all_pair_sync_scores_boot{kE}(insert_inds, insert_inds, kD, kB) = ...
                calc_pair_sync_scores(boot_trans_table, this_chans);
            
            % normalized mutual information
            [~, nmi_boot] = class_mut_info(horzcat(classes_shuffled{:}));
            nmi_boot(1:this_n_chans >= (1:this_n_chans)') = nan;
            mut_info_boot{kE}(insert_inds, insert_inds, kD, kB) = nmi_boot;
            
            % canonical correlation
            this_cca = nan(this_n_chans);
            for iC = 1:this_n_chans
                for jC = 1:iC-1
                    [~, ~, ps] = canoncorr(V_shuffled{iC}, V_shuffled{jC});
                    this_cca(iC, jC) = mean(ps);
                end
            end
            all_cca_boot{kE}(insert_inds, insert_inds, kD, kB) = this_cca;
        end        
        bootstrap_rstates{kE}{kD} = rstates;
    end
    
    % make transition sync CDF plot (comparison with bootstrap)
    cdf_trans = cum_trans_interp / cum_trans_interp(end);
    cdf_trans_boot = cum_trans_boot ./ cum_trans_boot(:, end);
    
    fh = figure;
    plot(cdf_interp_vals, cdf_trans);
    hold on;
    
    % Plot median and 95% CI of bootstrap CDFs with error bars
    boot_cdf_quantiles = quantile(cdf_trans_boot, [0.025, 0.5, 0.975]);
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
    
    % for matrices - find channels with no data
    nonempty_chan_combos = any(~isnan(mut_info_combined{kE}), 3);
    nonempty_cols = any(nonempty_chan_combos);
    nonempty_rows = any(nonempty_chan_combos');
    nonempty_chans = nonempty_cols | nonempty_rows;
    
    all_chan_names{kE} = all_chan_names{kE}(nonempty_chans);
    hr_chan_names{kE} = hr_chan_names{kE}(nonempty_chans);
    mut_info_combined{kE} = mut_info_combined{kE}(nonempty_chans, nonempty_chans, :);
    mut_info_boot{kE} = mut_info_boot{kE}(nonempty_chans, nonempty_chans, :, :);
    all_cca{kE} = all_cca{kE}(nonempty_chans, nonempty_chans, :);
    all_cca_boot{kE} = all_cca_boot{kE}(nonempty_chans, nonempty_chans, :, :);
    all_pair_sync_scores{kE} = all_pair_sync_scores{kE}(nonempty_chans, nonempty_chans, :);
    all_pair_sync_scores_boot{kE} = all_pair_sync_scores_boot{kE}(nonempty_chans, nonempty_chans, :, :);
end

%% Use bootstraps to compute z-scores
matrix_data = {mut_info_combined, all_cca, all_pair_sync_scores};
matrix_zscores = matrix_data;
boot_data = {mut_info_boot, all_cca_boot, all_pair_sync_scores_boot};

for kV = 1:length(matrix_data)
    for kE = 1:length(matrix_data{kV})
        for kD = 1:size(matrix_data{kV}{kE}, 3)
            mean_boot = mean(boot_data{kV}{kE}(:, :, kD, :), 4);
            std_boot = std(boot_data{kV}{kE}(:, :, kD, :), 0, 4);
            matrix_zscores{kV}{kE}(:, :, kD) = (matrix_data{kV}{kE}(:, :, kD) - mean_boot) ./ std_boot;
        end
    end
end

[mut_info_z, cca_z, pair_sync_scores_z] = matrix_zscores{:};

%%
save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'n_bootstrap', 'all_chan_names', 'hr_chan_names', ...
    'mut_info_combined', 'mut_info_boot', 'all_cca', 'all_cca_boot', ...
    'all_pair_sync_scores', 'all_pair_sync_scores_boot', 'bootstrap_rstates', ...
    'mut_info_z', 'cca_z', 'pair_sync_scores_z', '-v7.3');

%% Get channel permutation test distributions and p-values for contrasts
n_perm = 1e7;
contrasts = {{'SameRegion', 'CrossRegion'}, {'SameRegionNonL4', 'SameRegionL4'}};
matrix_zscores = {mut_info_z, cca_z, pair_sync_scores_z};
perm_pvals = zeros(length(matrix_zscores), length(contrasts), length(exp_info));

for kE = 1:length(exp_info)
    chan_names = all_chan_names{kE};
    n_chans = length(chan_names);

    perm_stats = zeros(length(matrix_zscores), length(contrasts), n_perm+1);
    for kP = 1:n_perm+1
        if kP == 1
            permutation = 1:n_chans;
        else
            permutation = randperm(n_chans);
        end
        
        chans_perm = chan_names(permutation);
        for kD = 1:length(matrix_zscores)
            z_perm_by_type = util.categorize_pair_data(matrix_zscores{kD}{kE}, chans_perm);
            for kC = 1:length(contrasts)
                perm_stats(kD, kC, kP) = mean(z_perm_by_type.(contrasts{kC}{1})) - ...
                    mean(z_perm_by_type.(contrasts{kC}{2}));
            end
        end
    end
    real_stats = perm_stats(:, :, 1);
    perm_stats = perm_stats(:, :, 2:end);
    
    % get conservative p_u which should be fine for our # of possible permutations
    % explained at https://arxiv.org/pdf/1603.05766.pdf
    perm_pvals(:, :, kE) = (sum(perm_stats > real_stats, 3) + 1) / (n_perm + 1);
end

%% Mutual information combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    
    % matrix plot
    med_mut_info = median(mut_info_combined{kE}, 3, 'omitnan');
    fh = plot_dist_mat(med_mut_info, hr_chan_names{kE}, ...
        sprintf('%s, median over %d days', this_info.type, this_ndays), ...
        'norm_mutual_info', [], all_chan_names{kE});
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_nmi_%s.fig', this_info.type)));
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_nmi_%s.svg', this_info.type)));
    
    % violin plot
    fh = make_violin(mut_info_combined{kE}, all_chan_names{kE}, {
        'SameRegion', 'CrossRegion', ...
        'SameRegionL4', 'SameRegionNonL4', ...
        'CrossRegionL4', 'CrossRegionNonL4'
        }, {
        'Within Region', 'Across Regions', ...
        {'Within Region,', 'Including L4'}, ...
        {'Within Region,', 'Excluding L4'}, ...
        {'Across Regions,', 'Including L4'}, ...
        {'Across Regions,', 'Excluding L4'}
        });
    fh.Position(3) = 950;
    
    % add dotted line separating plots
    hold on;
    myylims = get(gca, 'YLim');
    plot([2.5, 2.5], myylims, 'k--');
    
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
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('nmi_violin_%s.svg', this_info.type)));
end

%% CCA combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    % matrix plot
    med_cca = median(all_cca{kE}, 3, 'omitnan');    
    fh = plot_dist_mat(med_cca, hr_chan_names{kE}, ...
        sprintf('%s - median over %d days', this_info.type, this_ndays), ...
        'cca_mean_r', [], all_chan_names{kE});
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_cca_r_%s.fig', this_info.type)));
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_cca_r_%s.svg', this_info.type)));
    
    % violin plot
    fh = make_violin(all_cca{kE}, all_chan_names{kE}, {
        'SameRegion', 'CrossRegion', ...
        'SameRegionL4', 'SameRegionNonL4', ...
        'CrossRegionL4', 'CrossRegionNonL4'
        }, {
        'Within Region', 'Across Regions', ...
        {'Within Region,', 'Including L4'}, ...
        {'Within Region,', 'Excluding L4'}, ...
        {'Across Regions,', 'Including L4'}, ...
        {'Across Regions,', 'Excluding L4'}
        });
    fh.Position(3) = 950;
    
    % add dotted line separating plots
    hold on;
    myylims = get(gca, 'YLim');
    plot([2.5, 2.5], myylims, 'k--');
    
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
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('cca_violin_%s.svg', this_info.type)));
end

%% Transition pair synchrony combined plots

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    % matrix plot
    med_sync = median(all_pair_sync_scores{kE}, 3, 'omitnan');    
    fh = plot_dist_mat(med_sync, hr_chan_names{kE}, ...
        sprintf('%s - median over %d days', this_info.type, this_ndays), ...
        'trans_sync_scores', [], all_chan_names{kE});
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_trans_sync_%s.fig', this_info.type)));
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('med_trans_sync_%s.svg', this_info.type)));
    
    % violin plot
    fh = make_violin(all_pair_sync_scores{kE}, all_chan_names{kE}, {
        'SameRegion', 'CrossRegion', ...
        'SameRegionL4', 'SameRegionNonL4', ...
        'CrossRegionL4', 'CrossRegionNonL4'
        }, {
        'Within Region', 'Across Regions', ...
        {'Within Region,', 'Including L4'}, ...
        {'Within Region,', 'Excluding L4'}, ...
        {'Across Regions,', 'Including L4'}, ...
        {'Across Regions,', 'Excluding L4'}
        });
    fh.Position(3) = 950;

    % add dotted line separating plots
    hold on;
    myylims = get(gca, 'YLim');
    plot([2.5, 2.5], myylims, 'k--');
    
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
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', sprintf('trans_sync_violin_%s.svg', this_info.type)));
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
pair_sync_scores = nan(n_chans);

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

end

%% Functions to make violin plots

function [hf, vs] = make_violin(mats_over_days, chan_names, category_snames, category_hr_names)
% sname = struct name, hr_name = human-readable name (cells),
% use cell of char vecs or string array for multiline

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

vals_bytype = util.categorize_pair_data(mats_over_days, chan_names);
violin_s = struct;
my_colors = zeros(nplots, 3);
for kV = 1:nplots
    sname = category_snames{kV};
    violin_s.(sname) = vals_bytype.(sname);
    my_colors(kV, :) = violin_colors.(sname);
end

hf = figure;
vs = violinplot(violin_s);
% xticklabels(category_hr_names);

% use manual text labels to avoid tex markup and uneditable text
xticklabels([]);
ca = gca;
for k = 1:length(category_hr_names)
    labeltext = join(string(category_hr_names{k}), newline);
    text(k, ca.YLim(1), labeltext, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
end

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
real_vals_bytype = util.categorize_pair_data(real_mats, chan_names);
tval_real = tval_fn(real_vals_bytype);

% get permutation tvals
perm_vals_bytype = cellfun(@(perm_cca) util.categorize_pair_data(perm_cca, chan_names), ...
    num2cell(perm_mats, [1, 2, 3]));
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

