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
%             '2020-10-27'
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
%     struct('type', 'm1_v1_csd', 'all_chan_names', make_layer_names(["M1", "V1"]), ...
%         'days', {{
%             '2020-10-26_csd'
% %             '2020-10-27_csd'
%             '2020-10-28_csd'
%             '2020-10-29_csd'
%         }})
%     struct('type', 'bilateral_csd', 'all_chan_names', make_layer_names(["V1L", "V1R"]), ...
%         'days', {{
%             '2021-01-27_csd'
%             '2021-01-29_csd'
%             '2021-01-31_csd'
%             '2021-02-02_csd'
%         }})
    ];

exp_types = {exp_info.type}';
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

%% Get canonical correlation between score matrices 
score_cca(nmf_mfiles);

%% Compute pairwise NMI, canonical correlation, and transition synchrony, with shuffled versions

n_shuffle = 1000;
cdf_interp_vals = 0:0.01:1;

all_chan_names = cell(length(exp_info), 1); 
mut_info_combined = cell(length(exp_info), 1);
mut_info_shuffle = cell(length(exp_info), 1);
all_cca = cell(length(exp_info), 1);
all_cca_shuffle = cell(length(exp_info), 1);
all_pair_sync_scores = cell(length(exp_info), 1);
all_pair_sync_scores_shuffle = cell(length(exp_info), 1);

% see whether we are shuffling anew or using previous seeds
if ~exist('reuse_shuffle_seeds', 'var')
    reuse_shuffle_seeds = false;
end

if reuse_shuffle_seeds && ~exist('shuffle_seeds', 'var')
    warning('Flag to re-use shuffle seeds was set, but they are not in workspace');
    reuse_shuffle_seeds = false;
end

if ~reuse_shuffle_seeds
    shuffle_seeds = cell(length(exp_info), 1);
end

trans_sync_means_real = zeros(length(exp_info), 1);
trans_sync_means_shuffle = zeros(length(exp_info), n_shuffle);
trans_sync_pvals = zeros(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    all_chan_names{kE} = this_info.all_chan_names;
    n_chans = length(all_chan_names{kE});

    % gather distance information from each recording
    mut_info_combined{kE} = nan(n_chans, n_chans, this_ndays);
    mut_info_shuffle{kE} = nan(n_chans, n_chans, this_ndays, n_shuffle);
    all_cca{kE} = nan(n_chans, n_chans, this_ndays);
    all_cca_shuffle{kE} = nan(n_chans, n_chans, this_ndays, n_shuffle);
    all_pair_sync_scores{kE} = nan(n_chans, n_chans, this_ndays);
    all_pair_sync_scores_shuffle{kE} = nan(n_chans, n_chans, this_ndays, n_shuffle);
    
    if ~reuse_shuffle_seeds
        shuffle_seeds{kE} = cell(this_ndays, 1);
    end
    
    % for collecting interpolated cumulated frequency of transition sync scores across days
    cum_trans_interp = zeros(1, length(cdf_interp_vals));
    cum_trans_interp_shuffle = zeros(n_shuffle, length(cdf_interp_vals));
    
    % for computing transition sync means over days/shuffled runs
    total_trans_real = 0;
    total_trans_shuffle = zeros(1, n_shuffle);
    total_score_real = 0;
    total_score_shuffle = zeros(1, n_shuffle);
    
    for kD = 1:this_ndays
        this_day = this_info.input_s(kD).name;
        mt_paths = this_info.input_s(kD).mt_res_in;
        this_mfile = matfile(this_info.input_s(kD).nmf_res_out, 'Writable', true);
        this_chans = this_mfile.all_chans;
        this_n_chans = length(this_chans);
        
        % get "human readable channel names" i.e. depths
        mt_mfile = matfile(mt_paths{1});
        this_hr_chan_names = util.make_hr_chan_names(this_chans, mt_mfile.chan_locs);
                
        % compute transition sync scores
        trans_table = get_state_transitions(mt_paths, this_mfile);
        trans_table = calc_transition_synchronization(trans_table);
        
        % add to unnormalized CDF
        cum_edges = 0:1/this_n_chans:1;
        cum_vals = 0:1/(this_n_chans-1):1;
        cum_trans = histcounts(trans_table.sync_score, cum_edges, 'Normalization', 'cumcount');
        cum_trans_interp = cum_trans_interp + interp1(cum_vals, cum_trans, cdf_interp_vals);
        
        % for mean sync distribution
        total_trans_real = total_trans_real + height(trans_table);
        total_score_real = total_score_real + sum(trans_table.sync_score);
        
        % make individual transition plot
        fh = figure;
        plot_transitions(trans_table);
        title(sprintf('Transitions with sync scores - %s', this_day), 'Interpreter', 'none');
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('transitions_w_sync_%s.fig', this_day)));
        
        % also get pairwise mean scores
        pair_sync_scores = calc_pair_sync_scores(trans_table, this_chans);
        
        % make individual pairwise sync score plot
        fh = plot_dist_mat(pair_sync_scores, this_hr_chan_names, [], ...
            'trans_sync_scores', 'full_nodiag', this_chans);
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('pair_sync_scores_%s.fig', this_day)));
        % since this is symmetric and we don't want to duplicate samples,
        % nan-out the upper triangular part
        pair_sync_scores(1:this_n_chans >= (1:this_n_chans)') = nan;
        
        % compute mutual information
        classes = this_mfile.filtered_classes;
        classes = horzcat(classes{1}{:});
        [~, norm_mut_info] = class_mut_info(classes);
        
        % plot mutual information
        fh = plot_dist_mat(norm_mut_info, this_hr_chan_names, [], 'norm_mutual_info', ...
            'full_nodiag', this_chans);
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('norm_mut_info_%s.fig', this_day)));
        norm_mut_info(1:this_n_chans >= (1:this_n_chans)') = nan;
        
        % plot CCA
        cca_mat = this_mfile.all_cca_sim;
        cca_mat_sym = cca_mat;
        cca_mat_sym(isnan(cca_mat_sym)) = 0;
        cca_mat_sym = cca_mat_sym + cca_mat_sym';
        fh = plot_dist_mat(cca_mat_sym, this_hr_chan_names, [], 'cca', 'full_nodiag', this_chans);
        savefig(fh, fullfile(sr_dirs.results, this_day, sprintf('cca_%s_all.fig', this_day)));
        
        % save data to 3D array
        insert_inds = cellfun(@(c) find(strcmp(all_chan_names{kE}, c)), this_chans);
        mut_info_combined{kE}(insert_inds, insert_inds, kD) = norm_mut_info;
        all_cca{kE}(insert_inds, insert_inds, kD) = cca_mat;
        all_pair_sync_scores{kE}(insert_inds, insert_inds, kD) = pair_sync_scores;
        
        % also do bootstraps
        nmf_V = this_mfile.nmf_V;
        trans = this_mfile.filtered_transitions;
        models = cell(this_n_chans, 1);
        seeds = cell(n_shuffle, this_n_chans);

        for kS = 1:n_shuffle
            V_shuffled = cell(this_n_chans, 1);
            classes_shuffled = cell(this_n_chans, 1);

            for kC = 1:this_n_chans
                if reuse_shuffle_seeds
                    rng(shuffle_seeds{kE}{kD}{kS, kC});
                end
                
                [V_shuffled{kC}, classes_shuffled{kC}, models{kC}, seeds{kS, kC}] = ...
                    util.shuffle_scores_markov(nmf_V{1}{kC}, classes(:, kC), trans{1}{kC}, ...
                    models{kC}, ~reuse_shuffle_seeds);
            end
            
            % mean transition synchrony CDF
            % make struct to stand in for NMF mfile in get_state_transitions
            shuffle_nmf_info = struct;
            shuffle_nmf_info.nmf_classes = {classes_shuffled};
            shuffle_nmf_info.nmf_V = {V_shuffled};
            shuffle_nmf_info.time_axis = this_mfile.time_axis;
            shuffle_nmf_info.chan_names = this_mfile.chan_names;
            shuffle_trans_table = get_state_transitions(mt_paths, shuffle_nmf_info, ...
                struct('save_filtered_classes', false));
            shuffle_trans_table = calc_transition_synchronization(shuffle_trans_table);
            
            cum_trans = histcounts(shuffle_trans_table.sync_score, cum_edges, 'Normalization', 'cumcount');
            cum_trans_interp_shuffle(kS, :) = cum_trans_interp_shuffle(kS, :) + interp1(cum_vals, cum_trans, cdf_interp_vals);
            
            % for mean sync distribution
            total_trans_shuffle(kS) = total_trans_shuffle(kS) + height(shuffle_trans_table);
            total_score_shuffle(kS) = total_score_shuffle(kS) + sum(shuffle_trans_table.sync_score);
            
            % pair sync scores
            all_pair_sync_scores_shuffle{kE}(insert_inds, insert_inds, kD, kS) = ...
                calc_pair_sync_scores(shuffle_trans_table, this_chans);
            
            % normalized mutual information
            [~, nmi_shuffle] = class_mut_info(horzcat(classes_shuffled{:}));
            nmi_shuffle(1:this_n_chans >= (1:this_n_chans)') = nan;
            mut_info_shuffle{kE}(insert_inds, insert_inds, kD, kS) = nmi_shuffle;
            
            % canonical correlation
            this_cca = nan(this_n_chans);
            for iC = 1:this_n_chans
                for jC = 1:iC-1
                    [~, ~, ps] = canoncorr(V_shuffled{iC}, V_shuffled{jC});
                    this_cca(iC, jC) = mean(ps);
                end
            end
            all_cca_shuffle{kE}(insert_inds, insert_inds, kD, kS) = this_cca;
        end
        if ~reuse_shuffle_seeds
            shuffle_seeds{kE}{kD} = seeds;
        end
    end
    
    % make transition sync CDF plot (comparison with shuffled)
    cdf_trans = cum_trans_interp / cum_trans_interp(end);
    cdf_trans_shuffle = cum_trans_interp_shuffle ./ cum_trans_interp_shuffle(:, end);
    
    fh = figure;
    plot(cdf_interp_vals, cdf_trans);
    hold on;
    
    % Plot median and 95% CI of bootstrap CDFs with error bars
    shuffle_cdf_quantiles = quantile(cdf_trans_shuffle, [0.025, 0.5, 0.975]);
    xconf = [cdf_interp_vals, cdf_interp_vals(end:-1:1)];
    yconf = [shuffle_cdf_quantiles(3, :), shuffle_cdf_quantiles(1, end:-1:1)];
    plot(cdf_interp_vals, shuffle_cdf_quantiles(2, :), 'r');
    fill(xconf, yconf, 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    title({sprintf('Synchronization scores for all days vs. %d shuffled runs', n_shuffle), ...
        this_info.type}, 'Interpreter', 'none');
    xlabel('Sync score');
    ylabel('CDF');
    legend('Real', 'Shuffled', '95% CI', 'Location', 'northwest');
    figname = sprintf('transition_sync_cdf_%s', this_info.type);
    savefig(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.fig']));
    saveas(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.svg']));
    
    % mean sync score distribution
    trans_sync_means_real(kE) = total_score_real / total_trans_real;
    trans_sync_means_shuffle(kE, :) = total_score_shuffle ./ total_trans_shuffle;
    
    % pvalue based on normal fit (CLT)
    mean_mean = mean(trans_sync_means_shuffle(kE, :));
    std_mean = std(trans_sync_means_shuffle(kE, :));
    trans_sync_pvals(kE) = normcdf(trans_sync_means_real(kE), mean_mean, std_mean, 'upper');    
    
    % for matrices - find channels with no data
    nonempty_chan_combos = any(~isnan(mut_info_combined{kE}), 3);
    nonempty_cols = any(nonempty_chan_combos);
    nonempty_rows = any(nonempty_chan_combos');
    nonempty_chans = nonempty_cols | nonempty_rows;
    
    all_chan_names{kE} = all_chan_names{kE}(nonempty_chans);
    mut_info_combined{kE} = mut_info_combined{kE}(nonempty_chans, nonempty_chans, :);
    mut_info_shuffle{kE} = mut_info_shuffle{kE}(nonempty_chans, nonempty_chans, :, :);
    all_cca{kE} = all_cca{kE}(nonempty_chans, nonempty_chans, :);
    all_cca_shuffle{kE} = all_cca_shuffle{kE}(nonempty_chans, nonempty_chans, :, :);
    all_pair_sync_scores{kE} = all_pair_sync_scores{kE}(nonempty_chans, nonempty_chans, :);
    all_pair_sync_scores_shuffle{kE} = all_pair_sync_scores_shuffle{kE}(nonempty_chans, nonempty_chans, :, :);
end

save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'exp_types', 'n_shuffle', 'all_chan_names', ...
    'shuffle_seeds', 'exp_info', 'trans_sync_means_real', 'trans_sync_means_shuffle', 'trans_sync_pvals', '-v7.3'); 

%% Use bootstraps to compute z-scores

% make n_experiment x 1 struct array
pairwise_stats = struct(...
    'norm_mutual_info',     mut_info_combined, ...
    'cca',                  all_cca, ...
    'trans_sync_scores',    all_pair_sync_scores ...
    );

pairwise_stats_shuffle = struct(...
    'norm_mutual_info',     mut_info_shuffle, ...
    'cca',                  all_cca_shuffle, ...
    'trans_sync_scores',    all_pair_sync_scores_shuffle ...
    );


analysis_names = fieldnames(pairwise_stats_shuffle);

for kE = 1:length(exp_info)
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        mean_shuffle = mean(pairwise_stats_shuffle(kE).(aname), 4);
        std_shuffle = std(pairwise_stats_shuffle(kE).(aname), 0, 4);
        pairwise_stats(kE).([aname, '_z']) = (pairwise_stats(kE).(aname) - mean_shuffle) ./ std_shuffle;
        pairwise_stats(kE).([aname, '_pval']) = normcdf(pairwise_stats(kE).([aname, '_z']), 'upper');
        bonf_crit = 0.05 / sum(~isnan(pairwise_stats(kE).(aname)), [1, 2]);
        pairwise_stats(kE).([aname, '_sig_bonf']) = double(pairwise_stats(kE).([aname, '_pval']) < bonf_crit);
        pairwise_stats(kE).([aname, '_sig_bonf'])(isnan(pairwise_stats(kE).(aname))) = nan;
    end
end

save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'pairwise_stats', 'pairwise_stats_shuffle', ...
    'analysis_names', '-append');

%% Do (actual) bootstrap for bilateral vs. M1/V1 cross-region contrast
n_boot = 1e6;
boot_pvals_bilat_vs_m1v1 = struct;
ref_types = {'meansub'}; %, 'csd'};

for kR = 1:length(ref_types)
    rtype = ref_types{kR};
    res_mats{1} = pairwise_stats(kR*2 - 1); % M1/V1
    chan_names{1} = all_chan_names{kR*2 - 1};
    res_mats{2} = pairwise_stats(kR*2); % bilateral
    chan_names{2} = all_chan_names{kR*2};
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        cross_means = zeros(2, n_boot);
        actual_diff = 0;
        
        for kT = 1:2 % (type)
            n_days = size(res_mats{kT}.(aname), 3);
            
            actual_sum = 0;
            total_vals = 0;
            for kD = 1:n_days
                % trim matrix down to rows/cols present in this day
                mat = res_mats{kT}.(aname)(:, :, kD);
                empty_cols = all(isnan(mat));
                empty_rows = all(isnan(mat'));
                b_keep = ~empty_cols | ~empty_rows;
                this_chans = chan_names{kT}(b_keep);
                mat = mat(b_keep, b_keep);
                
                % find cross-region submatrix
                if kT == 1
                    chans1 = contains(this_chans, 'V1');
                    chans2 = contains(this_chans, 'M1');
                else
                    chans1 = contains(this_chans, 'V1R');
                    chans2 = contains(this_chans, 'V1L');
                end
                cross_mat = mat(chans1, chans2);
                actual_sum = actual_sum + sum(cross_mat, 'all');
                
                % iterate over bootstraps and accumulate sums of permuted matrix
                [m, n] = size(cross_mat);
                total_vals = total_vals + m*n;
                for kB = 1:n_boot
                    perm1 = randsample(m, m, true);
                    perm2 = randsample(n, n, true);
                    perm_mat = cross_mat(perm1, perm2);
%                     perm_mat = cross_mat(randsample(m*n, m*n, true)); % not legit - resampling whole matrix
                    cross_means(kT, kB) = cross_means(kT, kB) + sum(perm_mat, 'all');
                end
            end
            % convert sums to means
            cross_means(kT, :) = cross_means(kT, :) / total_vals;            
            if kT == 1
                actual_diff = actual_diff - actual_sum / total_vals;
            else
                actual_diff = actual_diff + actual_sum / total_vals;
            end
        end
        
        % collect bootstrap sample
        boot_diffs = diff(cross_means);
        
        % since this distribution is centered at the actual statistic value, want to find prob.
        % that it is > 2 * this value.
        boot_pvals_bilat_vs_m1v1.(rtype).(aname) = (sum(boot_diffs >= 2*actual_diff)+1)/(n_boot+1);
    end
end

save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'boot_pvals_bilat_vs_m1v1', '-append');

%% Get channel permutation test distributions and p-values for contrasts
n_perm = 1e7;

kC = 1;
contrasts(kC).name = 'SameVsCross';
contrasts(kC).type1 = 'SameRegion';
contrasts(kC).type2 = 'CrossRegion';
contrasts(kC).chans_to_shuffle = ""; % (all)

% region-specific
for reg = ["V1", "M1", "V1L", "V1R"]
    kC = kC + 1;
    contrasts(kC).name = sprintf('%sVsCross', reg);
    contrasts(kC).type1 = sprintf('In%s', reg);
    contrasts(kC).type2 = 'CrossRegion';
    contrasts(kC).chans_to_shuffle = "";
end

for reg = ["V1", "V1L", "V1R"]
    kC = kC + 1;
    contrasts(kC).name = sprintf('%s_NonL4VsL4', reg);
    contrasts(kC).type1 = sprintf('In%sNonL4', reg);
    contrasts(kC).type2 = sprintf('In%sL4', reg);
    contrasts(kC).chans_to_shuffle = reg + "_"; % to guarantee it's an exact match for V1
end

% Combined V1 L4 vs. non-L4 for bilateral
kC = kC + 1;
contrasts(kC).name = 'V1LR_NonL4VsL4';
contrasts(kC).type1 = 'InV1NonL4';
contrasts(kC).type2 = 'InV1L4';
contrasts(kC).chans_to_shuffle = "V1L/V1R";

shuffle_chansets = unique([contrasts.chans_to_shuffle], 'stable');
contrasts_by_shuffle_chanset = arrayfun(@(cs) contrasts([contrasts.chans_to_shuffle] == cs), ...
    shuffle_chansets, 'uni', false);

perm_test_res = cell2struct(cell(length(analysis_names) + 2, length(exp_info)), ...
    [analysis_names; {'seed'; 'exp_type'}], 1);

for kE = 1:length(exp_info)
    rng('shuffle');
    rstate = rng;
    perm_test_res(kE).seed = rstate.Seed;
    perm_test_res(kE).exp_type = exp_types{kE};
    
    % do days individually
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        trimmed_mats = cell(this_ndays, 1);
        trimmed_chans = cell(this_ndays, 1);
        for kD = this_ndays:-1:1
            % trim matrix down to rows/cols present in this day
            mat = pairwise_stats(kE).(aname)(:, :, kD);
            empty_cols = all(isnan(mat));
            empty_rows = all(isnan(mat'));
            b_keep = ~empty_cols | ~empty_rows;
            trimmed_chans{kD} = all_chan_names{kE}(b_keep);
            day_masks(kD) = util.get_symmetric_matrix_masks(trimmed_chans{kD});
            trimmed_mats{kD} = mat(b_keep, b_keep);
            % make symmetric for permutation convenience (and we don't care about the diagonal)
            n_kept = size(trimmed_mats{kD}, 1);
            trimmed_mats{kD}(1:n_kept > (1:n_kept)') = 0;
            trimmed_mats{kD} = trimmed_mats{kD} + trimmed_mats{kD}';
        end
        
        % loop through contrasts
        for kS = 1:length(shuffle_chansets)
            % split into independent permutation sets if necessary
            this_chansets = split(shuffle_chansets{kS}, "/");
            
            % skip contrasts relating to channels that aren't present
            if any(arrayfun(@(cs) ...
                    ~any(contains(all_chan_names{kE}, cs)), ...
                    this_chansets))
                continue;
            end
            
            n_chanset_contrasts = length(contrasts_by_shuffle_chanset{kS});
            b_keep_contrast = true(n_chanset_contrasts, 1);
            total_vals = zeros(n_chanset_contrasts, 2);
            for kC = 1:n_chanset_contrasts
                contrast = contrasts_by_shuffle_chanset{kS}(kC);
                % get real data and mean difference 
                vec1 = cell(this_ndays, 1);
                vec2 = cell(this_ndays, 1);
                
                for kD = 1:this_ndays
                    vec1{kD} = trimmed_mats{kD}(day_masks(kD).(contrast.type1));
                    vec2{kD} = trimmed_mats{kD}(day_masks(kD).(contrast.type2));
                end
                vec1 = cell2mat(vec1);
                vec2 = cell2mat(vec2);
                
                % skip contrasts relating to channel pair types that aren't present
                if isempty(vec1) || isempty(vec2)
                    b_keep_contrast(kC) = false;
                else
                    perm_test_res(kE).(aname).(contrast.name).data.pos = vec1;
                    perm_test_res(kE).(aname).(contrast.name).data.neg = vec2;
                    perm_test_res(kE).(aname).(contrast.name).real_diff = mean(vec1) - mean(vec2);
                    total_vals(kC, 1) = length(vec1);
                    total_vals(kC, 2) = length(vec2);
                end
            end
            n_chanset_contrasts = sum(b_keep_contrast);
            this_contrasts = contrasts_by_shuffle_chanset{kS}(b_keep_contrast);
            total_vals = total_vals(b_keep_contrast, :);
            
            % now do permutation test on relevant contrasts
            perm_contrast_diffs = zeros(n_chanset_contrasts, n_perm);
            
            for kD = 1:this_ndays
                % permute
                chans_to_permute = arrayfun(@(cs) find(contains(trimmed_chans{kD}, cs)), ...
                    this_chansets, 'uni', false);
                perm_lens = cellfun('length', chans_to_permute);
                permutation = 1:length(trimmed_chans{kD});
                
                for kP = 1:n_perm
                    for k = 1:length(chans_to_permute)
                        permutation(chans_to_permute{k}) = chans_to_permute{k}(randperm(perm_lens(k)));
                    end
                    permuted_mat = trimmed_mats{kD}(permutation, permutation);

                    for kC = 1:n_chanset_contrasts
                        contrast = this_contrasts(kC);
                        vec1 = permuted_mat(day_masks(kD).(contrast.type1));
                        vec2 = permuted_mat(day_masks(kD).(contrast.type2));
                        perm_contrast_diffs(kC, kP) = perm_contrast_diffs(kC, kP) + sum(vec1) / total_vals(kC, 1);
                        perm_contrast_diffs(kC, kP) = perm_contrast_diffs(kC, kP) - sum(vec2) / total_vals(kC, 2);
                    end
                end
            end
            
            % finally get pvalues
            for kC = 1:n_chanset_contrasts
                contrast = this_contrasts(kC);
                real_diff = perm_test_res(kE).(aname).(contrast.name).real_diff;
                pval = (sum(perm_contrast_diffs(kC, :) > real_diff) + 1) / (n_perm + 1);
                perm_test_res(kE).(aname).(contrast.name).pval = pval;
            end
        end
    end
end


save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'n_perm', 'perm_test_res', '-append');

%% Combined M1/V1 and bilateral V1 L4 vs. non-L4
combined_L4_perm_test_res = cell2struct(cell(length(analysis_names) + 2, length(exp_info)/2), ...
    [analysis_names; {'seed'; 'exp_type'}], 1);

for kExpPair = 1:length(exp_info)/2
    kEs = kExpPair * 2 + [-1, 0];
    this_infos = exp_info(kEs);
    this_ndayss = arrayfun(@(s) length(s.days), this_infos);
    
    rng('shuffle');
    rstate = rng;
    combined_L4_perm_test_res(kExpPair).seed = rstate.Seed;
    combined_L4_perm_test_res(kExpPair).exp_type = strjoin(exp_types(kEs), '-');

    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        trimmed_mats = cell(2, 1);
        trimmed_chans = cell(2, 1);
        day_masks = cell(2, 1);
        vec1 = cell(2, 1);
        vec2 = cell(2, 1);
        
        for kERel= 1:2
            this_ndays = this_ndayss(kERel);
            kE = kEs(kERel);
            trimmed_mats{kERel} = cell(this_ndays, 1);
            trimmed_chans{kERel} = cell(this_ndays, 1);
            vec1{kERel} = cell(this_ndays, 1);
            vec2{kERel} = cell(this_ndays, 1);
            
            for kD = this_ndays:-1:1
                % trim matrix down to rows/cols present in this day
                mat = pairwise_stats(kE).(aname)(:, :, kD);
                empty_cols = all(isnan(mat));
                empty_rows = all(isnan(mat'));
                b_keep = ~empty_cols | ~empty_rows;
                trimmed_chans{kERel}{kD} = all_chan_names{kE}(b_keep);
                day_masks{kERel}(kD) = util.get_symmetric_matrix_masks(trimmed_chans{kERel}{kD});
                trimmed_mats{kERel}{kD} = mat(b_keep, b_keep);

                % make symmetric for permutation convenience (and we don't care about the diagonal)
                n_kept = size(trimmed_mats{kERel}{kD}, 1);
                trimmed_mats{kERel}{kD}(1:n_kept > (1:n_kept)') = 0;
                trimmed_mats{kERel}{kD} = trimmed_mats{kERel}{kD} + trimmed_mats{kERel}{kD}';
                
                % grab real vectors of each type
                vec1{kERel}{kD} = trimmed_mats{kERel}{kD}(day_masks{kERel}(kD).InV1NonL4);
                vec2{kERel}{kD} = trimmed_mats{kERel}{kD}(day_masks{kERel}(kD).InV1L4);
            end
            vec1{kERel} = cell2mat(vec1{kERel});
            vec2{kERel} = cell2mat(vec2{kERel});
        end

        total_vals = zeros(1, 2); % number of non-L4 and L4 pairs, respectively

        vec1 = cell2mat(vec1);
        vec2 = cell2mat(vec2);
        total_vals(1) = length(vec1);
        total_vals(2) = length(vec2);

        combined_L4_perm_test_res(kExpPair).(aname).data.pos = vec1;
        combined_L4_perm_test_res(kExpPair).(aname).data.neg = vec2;
        real_diff = mean(vec1) - mean(vec2);
        combined_L4_perm_test_res(kExpPair).(aname).real_diff = real_diff;

        % now do permutation test
        perm_contrast_diffs = zeros(1, n_perm);
        
        for kERel= 1:2
            this_ndays = this_ndayss(kERel);
            for kD = 1:this_ndays
                % permute
                if kERel == 1
                    this_chansets = "V1";
                else
                    this_chansets = ["V1L", "V1R"];
                end
                
                chans_to_permute = arrayfun(@(cs) find(contains(trimmed_chans{kERel}{kD}, cs)), ...
                    this_chansets, 'uni', false);
                perm_lens = cellfun('length', chans_to_permute);
                permutation = 1:length(trimmed_chans{kERel}{kD});

                for kP = 1:n_perm
                    for k = 1:length(chans_to_permute)
                        permutation(chans_to_permute{k}) = chans_to_permute{k}(randperm(perm_lens(k)));
                    end
                    permuted_mat = trimmed_mats{kERel}{kD}(permutation, permutation);
                    vec1 = permuted_mat(day_masks{kERel}(kD).InV1NonL4);
                    vec2 = permuted_mat(day_masks{kERel}(kD).InV1L4);
                    perm_contrast_diffs(kP) = perm_contrast_diffs(kP) + sum(vec1) / total_vals(1);
                    perm_contrast_diffs(kP) = perm_contrast_diffs(kP) - sum(vec2) / total_vals(2);
                end
            end
        end

        % finally get pvalue
        pval = (sum(perm_contrast_diffs > real_diff) + 1) / (n_perm + 1);
        combined_L4_perm_test_res(kExpPair).(aname).pval = pval;
    end
end


save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'combined_L4_perm_test_res', '-append');


%% Some things for figures

analysis_names = fieldnames(pairwise_stats_shuffle);

hr_exp_names = struct(...
    'm1_v1',                    'M1/V1', ...
    'm1_v1_csd',                'M1/V1, CSD', ...
    'bilateral',                'Bilateral V1', ...
    'bilateral_csd',            'Bilateral V1, CSD', ...
    'm1_v1__bilateral',          'All recordings', ...
    'm1_v1_csd__bilateral_csd',  'All recordings, CSD');

hr_pair_type_names = struct(...
    'SameRegion',   'Within Region', ...
    'CrossRegion',  'Across Regions', ...
    'InV1L4',       {{'V1 L4 to', 'V1 non-L4'}}, ...
    'InV1NonL4',    {{'V1 non-L4 to', 'V1 non-L4'}});

analysis_titles = struct(...
    'norm_mutual_info', 'Normalized mutual information of classes', ...
    'cca',              'Canonical corrleation of NMF scores', ...
    'trans_sync_scores','Synchrony of class transitions');

analysis_ylabels = struct(...
    'norm_mutual_info', 'Normalized MI', ...
    'cca',              'Mean CCA r', ...
    'trans_sync_scores','Mean transition SYNC score');

violin_colors = struct(...
    'SameVsCross', [1, 0, 0; 0, 0, 1], ...
    'V1_NonL4VsL4', [219, 71, 28; 160, 32, 240] / 255); % 255, 0, 153]);    

%% Make violin plots, showing shuffled p-values for regional and L4 contrasts

% same vs. cross
type1 = 'SameRegion';
type2 = 'CrossRegion';

for kE = 1:length(exp_info)
    exp_name = exp_info(kE).type;
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};

        % build violin plot struct and other info
        contrast_s = perm_test_res(kE).(aname).SameVsCross;
        violin_s = struct;
        violin_s.(type1) = contrast_s.data.pos;
        violin_s.(type2) = contrast_s.data.neg;
        pvals = contrast_s.pval;
        
%         n_plots = 2;
%         my_colors = zeros(n_plots, 3);
%         xlabels = cell(n_plots, 1);
%         pvals = zeros(n_contrasts, 1);
%         violin_s = struct;
%         contrast_names = {contrasts.name};
% 
%         for kC = 1:n_contrasts
%             this_contrast_name = contrasts_to_use{kC};
%             my_colors(kC*2 + [-1, 0], :) = violin_colors.(this_contrast_name);
%             
%             this_contrast = contrasts(strcmp(contrast_names, this_contrast_name));
%             pos_type = this_contrast.type1;
%             neg_type = this_contrast.type2;
%             xlabels{kC*2 - 1} = hr_pair_type_names.(pos_type);
%             xlabels{kC*2} = hr_pair_type_names.(neg_type);
%             
%             contrast_s = perm_test_res(kE).(aname).(this_contrast_name);
%             violin_s.(pos_type) = contrast_s.data.pos;
%             violin_s.(neg_type) = contrast_s.data.neg;
%             pvals(kC) = contrast_s.pval;
%         end
        
        fh = figure;
        vs = violinplot(violin_s);
        
        % use manual text labels to avoid tex markup and uneditable text
        xticklabels([]);
        ca = gca;
        xlabels = {hr_pair_type_names.(type1), hr_pair_type_names.(type2)};
        for kV = 1:2
            labeltext = join(string(xlabels{kV}), newline);
            text(kV, ca.YLim(1), labeltext, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', 'FontName', 'Arial');    
            vs(kV).ViolinColor = violin_colors.SameVsCross(kV, :);
        end
        
        fh.Position = [300, 300, 400, 500];
        set(gca, 'FontName', 'Arial');

        title(sprintf('%s (%s)', analysis_titles.(aname), hr_exp_names.(exp_name)));
        ylabel(analysis_ylabels.(aname));
        
%         % add dotted line separating plots
%         hold on;
%         myylims = get(gca, 'YLim');
%         plot([2.5, 2.5], myylims, 'k--');
        
        % add p-value indicators
        pvals(pvals >= 0.05) = nan;
        hs = sigstar(mat2cell(1:length(pvals)*2, 1, repmat(2, 1, length(pvals))), pvals);
        set(hs(:, 2), 'VerticalAlignment', 'baseline', 'FontName', 'Arial', 'FontSize', 14);
        
        % make sure y axis goes high enough for sigstars
        ca.YLim(2) = max(ca.YLim(2), max(structfun(@max, violin_s)) * 1.2);
        
        figname = sprintf('violin_%s_%s', aname, exp_name);
        savefig(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.fig']));
        saveas(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.svg']));
    end
end

%% V1 L4 vs. non-L4, combined across M1/V1 and bilateral experiments
type1 = 'InV1NonL4';
type2 = 'InV1L4';

for kExpPair = 1:length(combined_L4_perm_test_res)
    exp_pair_name = combined_L4_perm_test_res(kExpPair).exp_type;
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        
        % build violin plot struct and other info
        contrast_s = combined_L4_perm_test_res(kExpPair).(aname);
        violin_s = struct;
        violin_s.(type1) = contrast_s.data.pos;
        violin_s.(type2) = contrast_s.data.neg;
        pvals = contrast_s.pval;
        
        fh = figure;
        vs = violinplot(violin_s);
        
        % use manual text labels to avoid tex markup and uneditable text
        xticklabels([]);
        ca = gca;
        xlabels = {hr_pair_type_names.(type1), hr_pair_type_names.(type2)};
        for kV = 1:2
            labeltext = join(string(xlabels{kV}), newline);
            text(kV, ca.YLim(1), labeltext, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', 'FontName', 'Arial');    
            vs(kV).ViolinColor = violin_colors.V1_NonL4VsL4(kV, :);
        end
        
        fh.Position = [300, 300, 400, 500];
        set(gca, 'FontName', 'Arial');
        
        title(sprintf('%s\n(combined)', analysis_titles.(aname)));
        ylabel(analysis_ylabels.(aname));
        
        % add p-value indicators
        pvals(pvals >= 0.05) = nan;
        hs = sigstar(mat2cell(1:length(pvals)*2, 1, repmat(2, 1, length(pvals))), pvals);
        set(hs(:, 2), 'VerticalAlignment', 'baseline', 'FontName', 'Arial', 'FontSize', 14);
        
        % make sure y axis goes high enough for sigstars
        ca.YLim(2) = max(ca.YLim(2), max(structfun(@max, violin_s)) * 1.2);
        
        figname = sprintf('violin_%s_%s', aname, exp_pair_name);
        savefig(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.fig']));
        saveas(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.svg']));
    end
end


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

% make it a full matrix
pair_sync_scores = pair_sync_scores + pair_sync_scores';

end

%% Function to shuffle labels on a pair of vectors to get permutation p-value
% function pval = perm_test(vec1, vec2, n_perm)
% % Tests the hypothesis that mean(vec1) > mean(vec2)
% 
% % coefficients for mean comparison:
% coefs = [ones(1, length(vec1)) / length(vec1), -ones(1, length(vec2)) / length(vec2)];
% vec_all = [vec1(:); vec2(:)];
% mean_diff_real = coefs * vec_all;
% mean_diff_perm = zeros(n_perm, 1);
% for kP = 1:n_perm
%     mean_diff_perm(kP) = coefs(randperm(length(coefs))) * vec_all;
% end
% 
% % See: https://arxiv.org/pdf/1603.05766.pdf
% pval = (sum(mean_diff_perm > mean_diff_real) + 1) / (n_perm + 1);
% 
% end

%% Function to make violin plots

% function [hf, vs] = make_violin(mats_over_days, chan_names, category_snames, category_hr_names)
% % sname = struct name, hr_name = human-readable name (cells),
% % use cell of char vecs or string array for multiline
% 
% violin_colors = struct(...
%     'SameChannel', [0.5 0.5 0.5], ...
%     'SameRegion', [0 0 1], ...
%     'SameRegionL4', [0.5 0.5 1], ...
%     'SameRegionNonL4', [0 0 0.5], ...
%     'CrossRegion', [1 0 0], ...
%     'CrossRegionL4', [1 0.5 0.5], ...
%     'CrossRegionNonL4', [0.5 0 0], ...
%     'InV1', [1, 0, 1], ...
%     'InV1L4', [1, 0.5, 1], ...
%     'InV1NonL4', [0.5, 0, 0.5], ...
%     'InM1', [0, 0.8, 0], ...
%     'InM1L4', [0.5, 1, 0.5], ...
%     'InM1NonL4', [0, 0.4, 0] ...
%     );
% 
% assert(length(category_snames) == length(category_hr_names), 'Mismatch in # of categories');
% nplots = length(category_snames);
% 
% vals_bytype = util.categorize_pair_data(mats_over_days, chan_names);
% violin_s = struct;
% my_colors = zeros(nplots, 3);
% for kV = 1:nplots
%     sname = category_snames{kV};
%     violin_s.(sname) = vals_bytype.(sname);
%     my_colors(kV, :) = violin_colors.(sname);
% end
% 
% hf = figure;
% vs = violinplot(violin_s);
% % xticklabels(category_hr_names);
% 
% % use manual text labels to avoid tex markup and uneditable text
% xticklabels([]);
% ca = gca;
% for k = 1:length(category_hr_names)
%     labeltext = join(string(category_hr_names{k}), newline);
%     text(k, ca.YLim(1), labeltext, 'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'top', 'FontName', 'Arial');
% end
% 
% for kV = 1:nplots
%     vs(kV).ViolinColor = my_colors(kV, :);
% end
% 
% end

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

