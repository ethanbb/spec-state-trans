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
hr_chan_names = cell(length(exp_info), 1); % human-readable
mut_info_combined = cell(length(exp_info), 1);
mut_info_shuffle = cell(length(exp_info), 1);
all_cca = cell(length(exp_info), 1);
all_cca_shuffle = cell(length(exp_info), 1);
all_pair_sync_scores = cell(length(exp_info), 1);
all_pair_sync_scores_shuffle = cell(length(exp_info), 1);
shuffle_seeds = cell(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);

    all_chan_names{kE} = this_info.all_chan_names;
    hr_chan_names{kE} = util.make_hr_chan_names(all_chan_names{kE}, 140);
    n_chans = length(all_chan_names{kE});

    % gather distance information from each recording
    mut_info_combined{kE} = nan(n_chans, n_chans, this_ndays);
    mut_info_shuffle{kE} = nan(n_chans, n_chans, this_ndays, n_shuffle);
    all_cca{kE} = nan(n_chans, n_chans, this_ndays);
    all_cca_shuffle{kE} = nan(n_chans, n_chans, this_ndays, n_shuffle);
    all_pair_sync_scores{kE} = nan(n_chans, n_chans, this_ndays);
    all_pair_sync_scores_shuffle{kE} = nan(n_chans, n_chans, this_ndays, n_shuffle);
    shuffle_seeds{kE} = cell(this_ndays, 1);
    
    % for collecting interpolated cumulated frequency of transition sync scores across days
    cum_trans_interp = zeros(1, length(cdf_interp_vals));
    cum_trans_interp_shuffle = zeros(n_shuffle, length(cdf_interp_vals));
    
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
        seeds = cell(n_shuffle, this_n_chans);

        for kS = 1:n_shuffle
            V_shuffled = cell(this_n_chans, 1);
            classes_shuffled = cell(this_n_chans, 1);

            for kC = 1:this_n_chans
                [V_shuffled{kC}, classes_shuffled{kC}, models{kC}, seeds{kS, kC}] = ...
                    util.shuffle_scores_markov(nmf_V{1}{kC}, classes(:, kC), trans{1}{kC}, models{kC});
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
        shuffle_seeds{kE}{kD} = seeds;
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
    
    % for matrices - find channels with no data
    nonempty_chan_combos = any(~isnan(mut_info_combined{kE}), 3);
    nonempty_cols = any(nonempty_chan_combos);
    nonempty_rows = any(nonempty_chan_combos');
    nonempty_chans = nonempty_cols | nonempty_rows;
    
    all_chan_names{kE} = all_chan_names{kE}(nonempty_chans);
    hr_chan_names{kE} = hr_chan_names{kE}(nonempty_chans);
    mut_info_combined{kE} = mut_info_combined{kE}(nonempty_chans, nonempty_chans, :);
    mut_info_shuffle{kE} = mut_info_shuffle{kE}(nonempty_chans, nonempty_chans, :, :);
    all_cca{kE} = all_cca{kE}(nonempty_chans, nonempty_chans, :);
    all_cca_shuffle{kE} = all_cca_shuffle{kE}(nonempty_chans, nonempty_chans, :, :);
    all_pair_sync_scores{kE} = all_pair_sync_scores{kE}(nonempty_chans, nonempty_chans, :);
    all_pair_sync_scores_shuffle{kE} = all_pair_sync_scores_shuffle{kE}(nonempty_chans, nonempty_chans, :, :);
end

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


analysis_names = fieldnames(pairwise_stats);

for kE = 1:length(exp_info)
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        mean_shuffle = mean(pairwise_stats_shuffle(kE).(aname), 4);
        std_shuffle = std(pairwise_stats_shuffle(kE).(aname), 0, 4);
        pairwise_stats(kE).([aname, '_z']) = (pairwise_stats(kE).(aname) - mean_shuffle) ./ std_shuffle;
    end
end

%%
save(fullfile(ssr_dirs.results, 'analysis_info.mat'), 'exp_types', 'n_shuffle', 'all_chan_names', ...
   'hr_chan_names', 'pairwise_stats', 'pairwise_stats_shuffle', 'shuffle_seeds', 'exp_info', '-v7.3');

%% Get channel permutation test distributions and p-values for contrasts
n_perm = 1e7;
contrasts(1).name = 'SameVsCross';
contrasts(1).type1 = 'SameRegion';
contrasts(1).type2 = 'CrossRegion';

contrasts(2).name = 'V1_NonL4VsL4';
contrasts(2).type1 = 'InV1NonL4';
contrasts(2).type2 = 'InV1L4';

perm_seeds = zeros(length(exp_info), 1);
perm_test_res = cell2struct(cell(length(analysis_names), length(exp_info)), analysis_names, 1);

% comparing bilateral to M1/V1
cross_exp_contrasts(1).name = 'Cross_BilatVsM1V1';
cross_exp_contrasts(1).type = 'CrossRegion';

cross_exp_perm_test_res = perm_test_res(1:2);
cross_exp_perm_test_res(1).exp_type = sprintf('%s_vs_%s', exp_types{[2,1]});
cross_exp_perm_test_res(2).exp_type = sprintf('%s_vs_%s', exp_types{[4,3]});

for kE = 1:length(exp_info)
    rng('shuffle');
    rstate = rng;
    perm_test_res(kE).seed = rstate.Seed;
    perm_test_res(kE).exp_type = exp_types{kE};
    
    chan_names = all_chan_names{kE};
    
    % for bilateral vs. m1/v1 (cross-experiment) contrasts
    ce_ind = floor((kE-1) / 2) + 1;
    if mod(kE, 2) == 0
        posneg = 'pos';
    else
        posneg = 'neg';
    end
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        stats_by_type = util.categorize_pair_data(pairwise_stats(kE).(aname), chan_names);
        for kC = 1:length(contrasts)
            cname = contrasts(kC).name;
            vec1 = stats_by_type.(contrasts(kC).type1);
            vec2 = stats_by_type.(contrasts(kC).type2);
            
            perm_test_res(kE).(aname).(cname).data.pos = vec1;
            perm_test_res(kE).(aname).(cname).data.neg = vec2;
            perm_test_res(kE).(aname).(cname).pval = perm_test(vec1, vec2, n_perm);
            
        end
        
        for kC = 1:length(cross_exp_contrasts)
            cname = cross_exp_contrasts(kC).name;
            vec = stats_by_type.(cross_exp_contrasts(kC).type);            
            cross_exp_perm_test_res(ce_ind).(aname).(cname).data.(posneg) = vec;
        end
    end
end

% get cross-exp p-values
for kE = 1:length(cross_exp_perm_test_res)
    rng('shuffle');
    rstate = rng;
    cross_exp_perm_test_res(kE).seed = rstate.Seed;
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};                
        for kC = 1:length(cross_exp_contrasts)
            cname = cross_exp_contrasts(kC).name;
            vec1 = cross_exp_perm_test_res(kE).(aname).(cname).data.pos;
            vec2 = cross_exp_perm_test_res(kE).(aname).(cname).data.neg;
            cross_exp_perm_test_res(kE).(aname).(cname).pval = perm_test(vec1, vec2, n_perm);
        end
    end
end

save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'n_perm', 'perm_test_res', 'cross_exp_perm_test_res', '-append');

%% Some things for figures

analysis_names = fieldnames(perm_pvals);

hr_exp_names = struct(...
    'm1_v1',        'M1/V1', ...
    'm1_v1_csd',    'M1/V1, CSD', ...
    'bilateral',    'Bilateral V1', ...
    'bilateral_csd','Bilateral V1, CSD');

analysis_titles = struct(...
    'norm_mutual_info', 'Normalized mutual information of classes', ...
    'cca',              'Canonical corrleation of NMF scores', ...
    'trans_sync_scores','Synchronization of class transitions');

analysis_ylabels = struct(...
    'norm_mutual_info', 'Normalized MI (z-score)', ...
    'cca',              'Mean CCA r (z-score)', ...
    'trans_sync_scores','Mean transition SYNC score (z-score)');

%% Make combined matrix plots (means of z-scores)

for kE = 1:length(exp_info)
    exp_name = exp_info(kE).type;
    exp_ndays = length(exp_info(kE).days);
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        zname = [aname, '_z'];
        mean_z = mean(pairwise_stats(kE).(zname), 3, 'omitnan');
        
        [fh, cb] = plot_dist_mat(mean_z, hr_chan_names{kE}, ...
            sprintf('%s, mean over %d days', hr_exp_names.(exp_name), exp_ndays), ...
            zname, 'lower_nodiag', all_chan_names{kE});
        
        title(sprintf('%s (%s)', analysis_titles.(aname), hr_exp_names.(exp_name)));
        cb.Label.String = 'z-score';
        cb.Label.FontSize = 11;
        
        figname = sprintf('mean_%s_%s', zname, exp_name);
        savefig(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.fig']));
        saveas(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.svg']));
    end
end

%% Make violin plots, showing shuffled p-values for regional and L4 contrasts
% for kE = 1:length(exp_info)
for kE = 1:2:length(exp_info)
    exp_name = exp_info(kE).type;
    exp_ndays = length(exp_info(kE).days);
    
    for kA = 1:length(analysis_names)
        aname = analysis_names{kA};
        zname = [aname, '_z'];
        fh = make_violin(pairwise_stats(kE).(zname), all_chan_names{kE}, ...
...             [contrasts.SameVsCross, contrasts.Same_L4VsNonL4], ...
...             {
...                 'Within Region'
...                 'Across Regions'
...                 {'Within Region', 'Outside L4'}
...                 {'Within Region', 'With L4'}
            [contrasts.SameVsCross, contrasts.V1_L4VsNonL4, contrasts.M1_L4VsNonL4], ...
            {
                'Within Region'
                'Across Regions'
                {'Within V1', 'Outside L4'}
                {'Within V1', 'With L4'}
                {'Within M1', 'Outside L4'}
                {'Within M1', 'With L4'}
            });
        fh.Position = [300, 300, 600, 500];
        set(gca, 'FontName', 'Arial');

        title(sprintf('%s (%s)', analysis_titles.(aname), hr_exp_names.(exp_name)));
        ylabel(analysis_ylabels.(aname));
        
        % add dotted line separating plots
        hold on;
        myylims = get(gca, 'YLim');
        plot([2.5, 2.5], myylims, 'k--');
        
        % add p-value indicators
        pval_s = perm_pvals(kE).(aname);
%         pvals = [pval_s.SameVsCross, pval_s.Same_L4VsNonL4];
        pvals = [pval_s.SameVsCross, pval_s.V1_L4VsNonL4, pval_s.M1_L4VsNonL4];
        pvals(pvals >= 0.05) = nan;
        hs = sigstar(mat2cell(1:length(pvals)*2, 1, repmat(2, 1, length(pvals))), pvals);
        set(hs(:, 2), 'VerticalAlignment', 'baseline', 'FontName', 'Arial', 'FontSize', 14);
        
%         figname = sprintf('violin_%s_%s', zname, exp_name);
%         savefig(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.fig']));
%         saveas(fh, fullfile(sr_dirs.results, 'res_figs', [figname, '.svg']));
    end
end

%% Finally, save
save(fullfile(sr_dirs.results, 'analysis_info.mat'), 'exp_info', '-append');


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

%% Function to shuffle labels on a pair of vectors to get permutation p-value
function pval = perm_test(vec1, vec2, n_perm)
% Tests the hypothesis that mean(vec1) > mean(vec2)

% coefficients for mean comparison:
coefs = [ones(1, length(vec1)) / length(vec1), -ones(1, length(vec2)) / length(vec2)];
vec_all = [vec1(:); vec2(:)];
mean_diff_real = coefs * vec_all;
mean_diff_perm = zeros(n_perm, 1);
for kP = 1:n_perm
    mean_diff_perm(kP) = coefs(randperm(length(coefs))) * vec_all;
end

% See: https://arxiv.org/pdf/1603.05766.pdf
pval = (sum(mean_diff_perm > mean_diff_real) + 1) / (n_perm + 1);

end

%% Function to make violin plots

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
    'CrossRegionNonL4', [0.5 0 0], ...
    'InV1', [1, 0, 1], ...
    'InV1L4', [1, 0.5, 1], ...
    'InV1NonL4', [0.5, 0, 0.5], ...
    'InM1', [0, 0.8, 0], ...
    'InM1L4', [0.5, 1, 0.5], ...
    'InM1NonL4', [0, 0.4, 0] ...
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
        'VerticalAlignment', 'top', 'FontName', 'Arial');
end

for kV = 1:nplots
    vs(kV).ViolinColor = my_colors(kV, :);
end

end


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

