% Analysis to decompose dynamics over all channels (after NMF) to look for
% components of global state, similar to Hudson et al.

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
            '2020-10-26'
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

% also get useful stuff from earlier analysis
info_s = load(fullfile(sr_dirs.results, 'analysis_info.mat'));

%% Do PCA on concatenated NMF scores

for kE = 1:length(exp_types)
    input_s = exp_info(kE).input_s;
    this_n_days = length(input_s);
    
    for kD = 1:this_n_days
        this_day_info = input_s(kD);
        nmf_mfile = matfile(this_day_info.nmf_res_out);
        
        scores = nmf_mfile.nmf_V;
        scores = scores{1};
        n_chans = length(scores);
        comps_per_chan = cellfun('size', scores, 2); % for mat2cell later
        full_dynamics = cell2mat(scores');
        fprintf('%s total NMF comps: %d\n', this_day_info.name, size(full_dynamics, 2));
        
        [pcs, pc_scores, eigvals] = pca(full_dynamics);
        
        % 2D histogram of data on 1st 2 PC axes
        fh = figure('Position', [200, 200, 530, 440]);
        ax_2d = subplot(4, 4, [2, 3, 4, 6, 7, 8, 10, 11, 12]);        
        [counts, xedges, yedges] = histcounts2(pc_scores(:, 1), pc_scores(:, 2), [100, 100]);
        xcenters = mean([xedges(1:end-1); xedges(2:end)]);
        ycenters = mean([yedges(1:end-1); yedges(2:end)]);
        sanePColor(xcenters, ycenters, imgaussfilt(counts', 3));
        pc1_lims = ax_2d.XLim;
        pc2_lims = ax_2d.YLim;
        
        % Marginal ksdensity plots
        ax_marginy = subplot(4, 4, [1, 5, 9]);
        ax_marginx = subplot(4, 4, [14, 15, 16]);
        
        [pc1_marg, xi] = ksdensity(pc_scores(:, 1));
        plot(ax_marginx, xi, pc1_marg, 'k');
        xlim(ax_marginx, pc1_lims);
        xlabel(ax_marginx, 'Prob. density');
        
        [pc2_marg, xi] = ksdensity(pc_scores(:, 2));
        plot(ax_marginy, pc2_marg, xi, 'k');
        ylim(ax_marginy, pc2_lims);
        ylabel(ax_marginy, 'Prob. density');
        
        ax_marginy.XDir = 'reverse';
        ax_marginx.YDir = 'reverse';
        ax_marginy.YAxis.Visible = false;
        ax_marginx.XAxis.Visible = false;
        ax_2d.FontName = 'Arial';
        ax_marginy.FontName = 'Arial';
        ax_marginx.FontName = 'Arial';
        
        savefig(fh, fullfile(sr_dirs.results, this_day_info.name, 'pc_space_data.fig'));
        
        % Plot PCA loadings in channel/frequency space
        loadings = nmf_mfile.nmf_U;
        loadings = loadings{1};
        freq_axis = nmf_mfile.freq_axis;
        loading_mat = zeros(length(freq_axis), n_chans);
        pcs_per_chan = mat2cell(pcs, comps_per_chan);
        pc_loadings_per_chan = cellfun(@(chan_loadings, chan_pcs) chan_loadings * chan_pcs, ...
            loadings, pcs_per_chan, 'uni', false);
        pc_loadings_all = cat(3, pc_loadings_per_chan{:}); % freq x pc# x chan
        pc_loadings_all = permute(pc_loadings_all, [1, 3, 2]); % freq x chan x pc#
        
        % get channel names
        mt_mfile = matfile(this_day_info.mt_res_in{1});
        this_chans = nmf_mfile.all_chans;
        this_hr_chan_names = util.make_hr_chan_names(this_chans, mt_mfile.chan_locs);
        
        % infer region labels (I know this is dumb)
        reg = split(this_chans, '_');
        reg = reg(:, 1);
        all_regs = unique(reg);
        k_reg = cellfun(@(r) find(strcmp(r, all_regs)), reg);
        
        n_to_plot = 6;
        fh = figure('Position', [11.5, 50, 1340, 730]);
        tl = tiledlayout(2, 3, 'TileSpacing', 'compact');
        
        for kP = 1:n_to_plot
            ax = nexttile;
            sanePColor(1:n_chans, freq_axis, pc_loadings_all(:, :, kP), false, true);
            ax.YScale = 'log';
            hold on;
            
            % Dividers between channels
            ylimits = get(gca, 'YLim');
            for kC = 2:n_chans
                edge_x = kC - 0.5;
                width = 1 + (k_reg(kC) ~= k_reg(kC-1));
                plot([edge_x, edge_x], ylimits, 'k', 'LineWidth', width);
            end
            xticks(1:n_chans);
            xticklabels(this_hr_chan_names);
            xtickangle(45);
            ylabel('Frequency (Hz)');
            title(sprintf('PC%d', kP));
            
            % put text for region names
            for kR = 1:length(all_regs)
                reg_center = mean(find(k_reg == kR));
                reg_name = all_regs{kR};
                text(reg_center, 0.03, reg_name, ...
                    'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold');
            end
            xlabel('Depth (um)');
        end
        title(tl, ['PCA loadings - ', this_day_info.name], 'Interpreter', 'none');
        savefig(fh, fullfile(sr_dirs.results, this_day_info.name, 'global_pca_loadings.fig'));
        
        % Variance explained - real and shuffled control
        
        frac_explained = cumsum(eigvals / sum(eigvals));
        n_vals_to_plot = max(8, min(sum(frac_explained < 0.95) + 1, length(frac_explained)));
        
        % Markov shuffle to compare to control
        n_shuffle = 1000;
        classes = nmf_mfile.nmf_classes;
        classes = classes{1};
        trans = nmf_mfile.nmf_transitions;
        trans = trans{1};
        frac_explained_shuffled = zeros(n_shuffle, length(frac_explained));

        scores_shuffled = cell(n_chans, 1);
        models = cell(n_chans, 1);
        for kS = 1:n_shuffle            
            for kC = 1:n_chans
                rng(info_s.shuffle_seeds{kE}{kD}{kS, kC});
                
                [scores_shuffled{kC}, ~, models{kC}] = util.shuffle_scores_markov( ...
                    scores{kC}, classes{kC}, trans{kC}, ...
                    models{kC}, false);
            end
            full_dynamics_shuffled = cell2mat(scores_shuffled');        
            [~, ~, eigvals_shuffled] = pca(full_dynamics_shuffled);
            frac_explained_shuffled(kS, :) = cumsum(eigvals_shuffled / sum(eigvals_shuffled));            
        end
        
        fh = figure;
        plot(frac_explained(1:n_vals_to_plot) * 100, 'b-s');
        hold on;
        
        % Plot median and 95% CI of shuffled frac_explained
        frac_explained_shuffled_quantiles = quantile(frac_explained_shuffled, [0.025, 0.5, 0.975]);
        xconf = [1:n_vals_to_plot, n_vals_to_plot:-1:1];
        yconf = [frac_explained_shuffled_quantiles(3, 1:n_vals_to_plot) * 100, ...
            frac_explained_shuffled_quantiles(1, n_vals_to_plot:-1:1) * 100];
        plot(frac_explained_shuffled_quantiles(2, 1:n_vals_to_plot) * 100, 'r-s');
        fill(xconf, yconf, 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        xlabel('Number of PCs');
        ylabel('Cumulative % variance');
        title(this_day_info.name, 'Interpreter', 'none');
        legend('Real', 'Shuffled', '95% CI', 'Location', 'northwest');
        savefig(fh, fullfile(sr_dirs.results, this_day_info.name, 'global_pca_explained.fig'));
    end
end
