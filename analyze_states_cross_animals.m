% Response to JNeurosci reviewer 2 question about consistency of states across recording
% sessions/animals. Two parts:
%  * cluster data combined from all animals, for a set of channels, and see how much each animal
%    contributes to each cluster (use same NMF method as for main text)
%  * plot some characteristic data from each global state (depending on result of first part)
%   - this is actually implemented in plot_example_data_from_classes.m.

% Channels to compare across animals
% silly thing to deal with the fact that the same region has a different name in different recordings
% channels = struct("M1V1", [
%         struct('tag', 'M1_L4', 'chan_name', 'M1_L4')
%         struct('tag', 'M1_Inf3', 'chan_name', 'M1_Inf3')
%         struct('tag', 'V1_L4', 'chan_name', 'V1_L4')
%         struct('tag', 'V1_Sup3', 'chan_name', 'V1_Sup3')
%     ], "BilatV1", [
%         struct('tag', 'V1_L4', 'chan_name', 'V1R_L4')
%         struct('tag', 'V1_Sup3', 'chan_name', 'V1R_Sup3')
%     ]);

rng('shuffle');

channels = struct(...
    "M1V1", struct('tag', 'V1_Sup2', 'chan_name', 'V1_Sup2'), ...
    "BilatV1", struct('tag', 'V1_Sup2', 'chan_name', 'V1R_Sup2'));

sr_dirs = prepSR;
cd(sr_dirs.script);
inds_filename = 'characteristic_inds.mat';

if ~exist(inds_filename, 'file')
    disp([inds_filename ' not found; running script to generate']);
    get_characteristic_class_inds;
end

inds_mfile = matfile(inds_filename);
days_by_type = struct("M1V1", categorical(inds_mfile.m1v1_dates), 'BilatV1', categorical(inds_mfile.bilatv1_dates));
alldays = inds_mfile.run_dates;
% topinds = inds_mfile.top_100_inds;
% pts_per_class = 100;

% Change to take a random sample of points instead
allinds = inds_mfile.characteristic_sortinds;
pts_per_class = 600;
n_reps = 500;

%% Load and combine data
% We end up with a big table

% do it in blocks to avoid running out of memory
reps_per_block = 50;
n_blocks = ceil(n_reps / reps_per_block);

mixed_state_day_nmi = cell(n_blocks, 1);
mixed_state_state_nmi = cell(n_blocks, 1);

for kB = 1:n_blocks
    rep_offset = (kB-1)*reps_per_block;
    this_n_reps = min(reps_per_block, n_reps - rep_offset);

    run_types = string(fieldnames(channels));
    null_classes = zeros(size(run_types));
    all_data = table(cell(length(run_types), 1), run_types, null_classes, ...
        'VariableNames', {'data', 'run_type', 'global_class'});

    for kT = 1:length(run_types)
        this_chans = channels.(run_types(kT));
        n_chans = length(this_chans);
        chan_ids = categorical({this_chans.tag});

        this_unique_days = days_by_type.(run_types(kT));
        this_n_days = length(this_unique_days);
        all_data.data{kT} = table(cell(this_n_days, 1), this_unique_days, 'VariableNames', {'data', 'day'});

        for kDay = 1:this_n_days
            this_day = char(this_unique_days(kDay));
            day_ind = find(strcmp(this_day, alldays));
    %         topinds_day = topinds{day_ind};
            allinds_day = allinds{day_ind};
            nmf_mfile = matfile(fullfile(this_day, 'nmf_res.mat'));
            freq_axis = nmf_mfile.freq_axis;
            chan_names = nmf_mfile.chan_names;
            chan_inds = cellfun(@(n) find(strcmp(n, chan_names)), {this_chans.chan_name});
            all_data.data{kT}.data{kDay} = table(cell(n_chans, 1), chan_ids(:), 'VariableNames', {'data', 'channel'});

            % load the actual data
            for kC = 1:n_chans
                pxx_cat_1chan = nmf_mfile.pxx_cat(chan_inds(kC), 1);
                pxx_cat_1chan = pxx_cat_1chan{1};

    %             inds_chan = topinds_day{kC};
                inds_chan = allinds_day{kC};
                classes = 1:length(inds_chan);
                classes(cellfun('length', inds_chan) < pts_per_class) = [];
                n_classes = length(classes);
                all_data.data{kT}.data{kDay}.data{kC} = table(cell(n_classes, 1), classes', 'VariableNames', {'data', 'class'});

                for kK = 1:n_classes
                    this_class = classes(kK);
                    this_inds = inds_chan{this_class};

                    ind_mat = arrayfun(@(k) this_inds(randsample(length(this_inds), pts_per_class)), 1:this_n_reps, 'uni', false);
                    data_mat3 = cell2mat(cellfun(@(inds) pxx_cat_1chan(:, inds)', reshape(ind_mat, 1, 1, []), 'uni', false));
                    all_data.data{kT}.data{kDay}.data{kC}.data{kK} = table(data_mat3, 'VariableNames', {'data'});

    %                 all_data.data{kT}.data{kDay}.data{kC}.data{kK} = table(pxx_cat_1chan(:, inds)', 'VariableNames', {'data'});
                end

                all_data.data{kT}.data{kDay}.data{kC} = explode_nested_tables(all_data.data{kT}.data{kDay}.data{kC});
            end
            all_data.data{kT}.data{kDay} = explode_nested_tables(all_data.data{kT}.data{kDay});
        end
        all_data.data{kT} = explode_nested_tables(all_data.data{kT});
    end
    all_data = explode_nested_tables(all_data);

%     break;
    
    % Now for each channel of interest, do NMF on whole concatenated dataset and compute
    % mean nmi between mixed states and both local states and animals.

    % res_mfile = matfile('cross_animal_states.mat', 'Writable', true);
    % xval_fig_dir = fullfile('res_figs', 'cross_animal_nmf_xval');
    cois = categories(all_data.channel);
    n_coi = length(cois);
    % n_comps = nan(n_coi, 1);

    % break apart by channel
%     each_chan_data = cellfun(@(c) all_data(all_data.channel == c, :), cois, 'uni', false);

    % Us = cell(n_coi, 1);
    
    this_mixed_state_day_nmi = nan(this_n_reps, n_coi);
    this_mixed_state_state_nmi = nan(this_n_reps, n_coi);

    for kC = 1:n_coi
        this_chan_data = all_data(all_data.channel == cois{kC}, :);
        this_unique_days = unique(this_chan_data.day);
        this_local_classes = arrayfun(@(d) this_chan_data.class(this_chan_data.day == d), ...
            this_unique_days, 'uni', false);
        
        this_chan_spectra = this_chan_data.data;
        this_chan_days = this_chan_data.day;
        
        classes_example = {};
        parfor kR = 1:this_n_reps
            data = this_chan_spectra(:, :, kR);
            
            this_n_comps = util.nmf_ncomps_xval(data, 1, 15, 0.01);

    %         figure(hfig);
    %         title(sprintf('NMF cross-validation (%s)', cois{kC}));
    %         savefig(hfig, fullfile(xval_fig_dir, ['nmf_xval_', cois{kC}, '.fig']));
    %         n_comps(kC) = this_n_comps;
    %         close(hfig);

            % do NMF with inferred # of components
            [V, U] = sp_nnmf(data, this_n_comps, [], [], 500000);

%             % do sorting and normalization as in concat_and_nmf (should probably be encapuslated)
%             [~, peak_freqinds] = max(U);
%             [~, order] = sort(peak_freqinds);
%             U = U(:, order);
%             V = V(:, order);

            % normalize
            norm_factor = vecnorm(U);
            U = U ./ norm_factor;
            V = V .* norm_factor;

            % get most likely "class"
            [~, classes] = max(V, [], 2);
            
            if kR == 1
                classes_example{kR} = classes;
            end
            
            % day NMI
            [~, this_mixed_state_day_nmi(kR, kC)] = class_mut_info(classes, this_chan_days);
            
            % mean individual vs. mixed class nmi
            [~, local_vs_mixed_nmi] = arrayfun(@(kD) ...
                class_mut_info(classes(this_chan_days == this_unique_days(kD)), this_local_classes{kD}), ...
                1:this_n_days);
            this_mixed_state_state_nmi(kR, kC) = mean(local_vs_mixed_nmi);

    %         Us{kC} = U;
    %         this_data.nmf_V = V;
    %         this_data.global_class = classes;
    %         each_chan_data{kC} = this_data;
        end
        all_data.global_class(all_data.channel == cois{kC}) = classes_example{1};
    end
    
    mixed_state_day_nmi{kB} = this_mixed_state_day_nmi;
    mixed_state_state_nmi{kB} = this_mixed_state_state_nmi;
end

mixed_state_day_nmi = cell2mat(mixed_state_day_nmi);
mixed_state_state_nmi = cell2mat(mixed_state_state_nmi);

% also save table with 1 example repetition
single_rep_data = all_data;
single_rep_data.data = single_rep_data.data(:, :, 1);

save('mixed_state_nmi.mat', 'mixed_state_day_nmi', 'mixed_state_state_nmi', 'single_rep_data', '-v7.3');
 
%% Stick plot comparing NMI

load('mixed_state_nmi.mat');

% violin_s.day_to_mixed_state = mixed_state_day_nmi;
% violin_s.individual_state_to_mixed_state = mixed_state_state_nmi;

data_quantiles = zeros(2, 3);
quantile_ps = [0.025, 0.5, 0.975];

data_quantiles(1, :) = quantile(mixed_state_day_nmi, quantile_ps);
data_quantiles(2, :) = quantile(mixed_state_state_nmi, quantile_ps);

neg = data_quantiles(:, 2) - data_quantiles(:, 1);
pos = data_quantiles(:, 3) - data_quantiles(:, 2);

diff_dist = mixed_state_state_nmi - mixed_state_day_nmi;
mean_diff = mean(diff_dist);
diff_dist_0centered = diff_dist - mean_diff;
pval = (sum(diff_dist_0centered >= mean_diff) + 1) / (length(diff_dist) + 1);
pval(pval > 0.05) = nan;

figure;
% vs = violinplot(violin_s);
% vs(2).ViolinColor = 'r';
errorbar([1, 2], data_quantiles(:, 2), neg, pos, 'k.', 'LineStyle', 'none', ...
    'MarkerSize', 10, 'LineWidth', 1, 'CapSize', 20);
xlim([0.7, 2.3]);
ylim([0, 0.3]);

hs = sigstar({1:2}, pval);
set(hs(:, 2), 'VerticalAlignment', 'baseline', 'FontName', 'Arial', 'FontSize', 14);

xticks([1, 2]);
xticklabels({'Animal ID', 'Per-animal state'});
ylabel('Normalized mutual information');
title('Information conveyed by mixed states');

%%
% res_mfile.chan_names = cois;
% res_mfile.n_comps = n_comps;
% res_mfile.nmf_Us = Us;
% res_mfile.data_by_chan = each_chan_data;
% res_mfile.freq_axis = freq_axis;

%% For each channel make pie charts of days included in each global class.

% cois = res_mfile.chan_names;
% n_comps = res_mfile.n_comps;
% Us = res_mfile.nmf_Us;
% each_chan_data = res_mfile.data_by_chan;
% freq_axis = res_mfile.freq_axis;

% reuse a single one of the random samples
cois = categories(single_rep_data.channel);
n_chans = length(cois);

% figure;
% h_scatter = axes;
% h_scatter.YScale = 'Log';
% xlabel('Relative entropy of day distribution');
% ylabel('Peak global class frequency (Hz)');
% title('Entropies of global classes vs. frequency');
% xlim([0, 1]);
% ylim([freq_axis(1), freq_axis(end)]);
% grid on;
% hold on;
% legend('Interpreter', 'none');

for kC = 1:n_chans
%     chan_data = each_chan_data{kC};
    chan_data = single_rep_data(single_rep_data.channel == cois{kC}, :);
    chan_classes = unique(chan_data.global_class);
    n_classes = length(chan_classes);
    chan_day_data = removecats(chan_data.day);
    days = unique(chan_day_data);
    this_n_days = length(days);
%     U = Us{kC};
    
    figure;
    n_cols = 3;
    n_rows = ceil(n_classes / 3);
    
    
%     tl = tiledlayout(n_cols, n_rows, 'Padding', 'compact', 'TileSpacing', 'compact');
%     title(tl, sprintf('Rats contributing to example\ncross-animal mixed states'));
    
    % n per class, for resizing pies
    n_per_class = sum(chan_data.global_class == chan_classes(:)');
    n_per_class_rel = n_per_class ./ max(n_per_class);
    scaled_area_side = sqrt(n_per_class_rel);
    
    for kK = 1:n_classes
        this_class = chan_classes(kK);
        class_days = chan_day_data(chan_data.global_class == this_class);
%         nexttile;
        subplot(n_cols, n_rows, kK);
        p = pie(class_days);
        
        pText = findobj(p, 'Type', 'text');
        set(pText, 'Visible', 'off');  % turn off labels
        pWedge = findobj(p, 'Type', 'patch');
        for kW = 1:length(pWedge)
            set(pWedge(kW), 'Vertices', get(pWedge(kW), 'Vertices') * scaled_area_side(kK));
        end
        
%         title(sprintf('Mixed state %d', kK));
    end
    
    suptitle(sprintf('Rats contributing to example\ncross-animal mixed states'));
    
%     % Stacked bar plots of individual channel states in each mixed state
%     figure;
%     tl2 = tiledlayout(n_cols, n_rows, 'Padding', 'compact', 'TileSpacing', 'compact');
%     title (tl2, sprintf('Rats and rat-specific states contributing to global classes (%s)', cois{kC}), 'Interpreter', 'none');
%     max_local_classes = max(chan_data.class);
%     
%     for kK = 1:n_classes
%         this_class = chan_classes(kK);
%         day_class_mat = zeros(this_n_days, max_local_classes);
%         for kD = 1:this_n_days
%             local_classes_inday = chan_data.class(chan_day_data == days(kD) & chan_data.global_class == this_class);
%             day_class_mat(kD, :) = histcounts(local_classes_inday, 0.5:max_local_classes+0.5);
%         end
%         nexttile;
%         bar(day_class_mat, 'stacked');
%         xticks(1:this_n_days);
%         xlabel('Rat #');
%         ylabel('# of windows');
%         title(sprintf('Mixed state %d', kK));
%     end
    
%     % Plot global states loadings
%     figure;
%     sanePColor(1:n_classes, freq_axis, U(:, chan_classes), false, true);
%     set(gca, 'YScale', 'log');
%     title(sprintf('Loadings for %s', cois{kC}), 'Interpreter', 'none');
%     xlabel('Global class');
%     ylabel('Frequency (Hz)');
    
%     % Alternative - plot as matrix of probabilities/distributions?
%     day_dist_mat = zeros(this_n_days, n_classes);
%     for kDay = 1:this_n_days
%         day_classes = chan_data.global_class(chan_day_data == days(kDay));
%         for kK = 1:n_classes
%             day_dist_mat(kDay, kK) = sum(day_classes == chan_classes(kK));
%         end
%     end
%     day_dist_mat = day_dist_mat ./ sum(day_dist_mat);
    
%     figure;
%     sanePColor(day_dist_mat);
%     colormap('jet');
%     xlabel('Global class');
%     ylabel('Frac. in each day');
%     title(sprintf('Global class distributions - %s', cois{kC}), 'Interpreter', 'none');
    
%     % plot entropy by peak frequency of each class
%     [~, peak_freqinds] = max(U(:, chan_classes));
%     peak_freqs = freq_axis(peak_freqinds);
%     entropies = -sum(day_dist_mat .* log2(day_dist_mat + eps));
%     max_entropy = log2(this_n_days);
%     rel_entropies = entropies / max_entropy;
%     scatter(h_scatter, rel_entropies, peak_freqs, 'filled', 'DisplayName', cois{kC});

%     % NMI between local and global classes
%     global_nmis = zeros(this_n_days, 1);
%     for kDay = 1:this_n_days
%         day_global_classes = chan_data.global_class(chan_day_data == days(kDay));
%         day_local_classes = chan_data.class(chan_day_data == days(kDay));
%         [~, global_nmis(kDay)] = class_mut_info(day_global_classes, day_local_classes);
%     end
%     
%     figure;
%     boxplot(global_nmis);
%     xticks([]);
%     ylabel('NMI between per-animal and cross-animal states');
end

