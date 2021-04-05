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
    med_mut_info = median(mut_info_combined, 3, 'omitnan');
    
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

chan_names = cell(length(exp_info), 1);
all_dists = cell(length(exp_info), 1);
all_dists_null = cell(length(exp_info), 1);

for kE = 1:length(exp_info)
    this_info = exp_info(kE);
    this_ndays = length(this_info.days);
    nmf_res_paths = {this_info.input_s.nmf_res_out}';
    nmf_res_mfiles = cellfun(@matfile, nmf_res_paths, 'uni', false);

    chan_names{kE} = this_info.all_chan_names;
    n_chans = length(chan_names{kE});

    % gather distance information from each recording
    all_dists{kE} = nan(n_chans, n_chans, this_ndays);
    all_dists_null{kE} = nan(n_chans, n_chans, this_ndays);

    for kD = 1:this_ndays
        insert_inds = cellfun(@(c) find(strcmp(chan_names{kE}, c)), nmf_res_mfiles{kD}.chan_names);

        this_mfile = nmf_res_mfiles{kD};
        all_dists{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dists;
        all_dists_null{kE}(insert_inds, insert_inds, kD) = this_mfile.allL2_dists_null;
    end

    med_dists = median(all_dists{kE}, 3, 'omitnan');

    % eliminate channels with no data
    chans_empty = all(isnan(med_dists)) & all(isnan(med_dists'));
    chan_names{kE} = chan_names{kE}(~chans_empty);
    med_dists = med_dists(~chans_empty, ~chans_empty);
    all_dists{kE} = all_dists{kE}(~chans_empty, ~chans_empty, :);
    all_dists_null{kE} = all_dists_null{kE}(~chans_empty, ~chans_empty, :);

    % make human-readable channel names
    hr_chan_names = util.make_hr_chan_names(chan_names{kE}, 140);

    hf = plot_dist_mat(med_dists, hr_chan_names, ...
        sprintf('%s - median over %d days', this_info.type, this_ndays), 'L2_dist');
    savefig(hf, fullfile(sr_dirs.results, 'res_figs', sprintf('med_score_dist_%s.fig', this_info.type)));
    
    % Same thing but with difference from null model
    med_dists_from_null = median(all_dists_null{kE} - all_dists{kE}, 3, 'omitnan');
    hf = plot_dist_mat(med_dists_from_null, hr_chan_names, ...
        sprintf('%s (median)', this_info.type), 'L2_dist_from_null');
    savefig(hf, fullfile(sr_dirs.results, 'res_figs', sprintf('med_score_dist_fromnull_%s.fig', this_info.type)));
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

%% Make violin plot

types = {'Single channel replication', 'Within region', 'Across regions'};
type_snames = cellfun(@matlab.lang.makeValidName, types, 'uni', false);
n_types = length(types);

for kE = 1:length(exp_info)
    this_chan_names = chan_names{kE};
    n_chans = length(this_chan_names);
    this_regions = strtok(this_chan_names, '_');
    same_region = this_regions == this_regions';
    
    % get indices of each type of pair
    linear_inds = cell(n_types, 1);
    linear_inds{1} = find(eye(n_chans));
    linear_inds{2} = setdiff(find(same_region), linear_inds{1});
    linear_inds{3} = find(~same_region);

    dists_bytype = struct;
    dists_null_bytype = struct;
    dists_frac_of_null = struct;
    
    for kT = 1:n_types
        %mean_kldiv_oftype = @(mat) mean(mat(linear_inds{kT}));
        tname = type_snames{kT};

        dists_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_dists{kE}, [1, 2]), 'uni', false);
        dists_bytype.(tname) = vertcat(dists_perday{:});

        dists_null_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_dists_null{kE}, [1, 2]), 'uni', false);
        dists_null_bytype.(tname) = vertcat(dists_null_perday{:});

        dists_frac_of_null.(tname) = dists_bytype.(tname) ./ dists_null_bytype.(tname);
    end

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

    hf = figure;
    viol = violinplot(dists_frac_of_null);
    xticklabels(types);
    ylabel('Score distance (frac. of null model)');
    title({'Euclidean distance of aligned NMF scores, compared to null model', ...
        exp_info(kE).type}, 'Interpreter', 'none');
    savefig(hf, fullfile(sr_dirs.results, 'res_figs', ...
        sprintf('dist_vs_null_violin_%s.fig', exp_info(kE).type)));

    [h, p] = ttest2(dists_frac_of_null.(type_snames{3}), dists_frac_of_null.(type_snames{2}), ...
        0.01, 'right');

    fprintf('%s: p = %g\n', exp_info(kE).type, p);
end
