% Cluster states for each channel on each day using non-negative matrix factorization, then
% generate null model data, find KL divergence between pairs of channels (after transformation),
% and make plots.

sr_dirs = prepSR;

days = {
    '2020-01-30'
    '2020-01-31'
    '2020-02-06'
    '2020-03-05'
    '2020-03-06'
    '2020-03-10'
    '2020-03-11'
    };

n_days = length(days);

%% Collect datasets on each day
input_s = struct('name', days, ...
    'mt_res_in', cell(n_days, 1), ...
    'nmf_res_out', cell(n_days, 1), ...
    'xval_fig_dir', cell(n_days, 1));

for kD = 1:n_days
    % Build input to concat_and_nmf
    
    curr_day = days{kD};
    input_s(kD).nmf_res_out = fullfile(sr_dirs.results, curr_day, 'nmf_res.mat');
    input_s(kD).xval_fig_dir = fullfile(sr_dirs.results, curr_day);
    
    % find all "layers" results files under this day
    res_fn = 'mt_res_layers.mat';
    res_entries = dir(fullfile(sr_dirs.results, curr_day, '*', res_fn));
    dirs = sort({res_entries.folder});
    input_s(kD).mt_res_in = fullfile(dirs, res_fn);
end

%% Do NMF
nmf_mfiles = concat_and_nmf(input_s);

%% Make null model data
gen_null_model_data(nmf_mfiles);

%% Do KL divergence analysis
kl_divergence_analysis(nmf_mfiles);


%% Get and plot median KL divergence for each channel pair (here still assuming V1 super...etc. labels are comparable)
% For each day, use only chan names specified in the csd results file (since I went back and
% filtered out channels that don't actually correspond to the layer they're supposed to be)

layers = {'L2/3', 'L4', 'L5'};
chans = [strcat('V1_', layers), strcat('M1_', layers)]; % how they will be plotted
n_chans = length(chans);

all_mean_kl_divs = nan(n_chans, n_chans, n_days);
all_mean_kl_divs_null = zeros(n_chans, n_chans, n_days);

for kD = 1:n_days
    res_kld = load(fullfile(sr_dirs.results, days{kD}, 'nmf_res.mat'), 'kl_divs', 'kl_divs_null', 'chan_names');
    csd_chans_V1_s = load(fullfile(sr_dirs.results, days{kD}, 'csd_V1.mat'), 'chan_names');
    csd_chans_MC_s = load(fullfile(sr_dirs.results, days{kD}, 'csd_MC.mat'), 'chan_names');
    
    data_cnames = res_kld.chan_names;
    
    good_chan_set = [strcat('V1_', csd_chans_V1_s.chan_names), strcat('M1_', csd_chans_MC_s.chan_names)];
    b_take = ismember(data_cnames, good_chan_set);
    data_good_cnames = data_cnames(b_take);
    insert_inds = cellfun(@(c) find(strcmp(chans, c)), data_good_cnames);
    
    all_mean_kl_divs(insert_inds, insert_inds, kD) = mean(res_kld.kl_divs(b_take, b_take, :), 3);
    all_mean_kl_divs_null(insert_inds, insert_inds, kD) = mean(res_kld.kl_divs_null(b_take, b_take, :), 3);
end

med_kl_divs = nanmedian(all_mean_kl_divs, 3);
hf = plot_kldiv_mat(med_kl_divs, chans, sprintf('Median over %d days', n_days));
savefig(hf, fullfile(sr_dirs.results, 'res_figs', 'med_kl_div_days.fig'));

%% Try same thing but with difference from null model

med_kl_div_from_null = nanmedian(all_mean_kl_divs_null - all_mean_kl_divs, 3);
hf = plot_kldiv_mat(med_kl_div_from_null, chans, 'Amount below divergence from null model (median)');
savefig(hf, fullfile(sr_dirs.results, 'res_figs', 'med_kl_div_fromnull_days.fig'));

% make a graph plot out of it
g = digraph(med_kl_div_from_null);
lwidths = 5 * g.Edges.Weight / max(g.Edges.Weight);
edge_labels = arrayfun(@(w) sprintf('%.2f', w), g.Edges.Weight, 'uni', false);

hf_graph = figure;
plot(g, 'LineWidth', lwidths_frac, 'NodeLabel', chans, 'EdgeLabel', edge_labels, ...
     'Interpreter', 'none', 'NodeFontSize', 12, 'NodeFontWeight', 'bold', 'Layout', 'layered', ...
     'Direction', 'right', 'Sources', 1:3, 'Sinks', 4:6)

title(sprintf('Synchrony of channel pairs\n(bits below null model KL divergence)'));
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
    
    kl_divs_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_mean_kl_divs, [1, 2]), 'uni', false);
    kl_divs_bytype.(type_snames{kT}) = vertcat(kl_divs_perday{:});
    
    kl_divs_null_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_mean_kl_divs_null, [1, 2]), 'uni', false);
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
