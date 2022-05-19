% Generate panels for figures 3 and 4.

sr_dirs = prepSR;

sample_day = '2020-10-26';
sample_time = '13-31-00';
sample_region = 'M1';
switch sample_region
    case {'M1', 'V1L'}
        sample_probe = 'Probe1';
        chan_offset = 0;
    case {'V1', 'V1R'}
        sample_probe = 'Probe2';
        chan_offset = 64;
end
sample_channames = {'Inf3', 'Inf4'};
sample_fullchannames = strcat(sample_region, '_', sample_channames);
sample_timerange = [11400, 11700];
figpos = [0, 0, 784, 238];

mt_res_first = matfile(fullfile(sr_dirs.results, sample_day, '11-30-00', 'mt_res_layers.mat'));
mt_options = mt_res_first.options;
winlen = mt_options.window;
time_offset = mt_res_first.time_grid(1, end) + winlen / 2;

mt_res = matfile(fullfile(sr_dirs.results, sample_day, sample_time, 'mt_res_layers.mat'));
timerange_offset = sample_timerange - time_offset;
time_mask = mt_res.time_grid >= timerange_offset(1) & mt_res.time_grid <= timerange_offset(2);

xaxis = linspace(0, sample_timerange(2) - sample_timerange(1), sum(time_mask));

% Channels
mfile_csd = matfile(fullfile(sr_dirs.results, sample_day, sprintf('csd_%s.mat', sample_region)));
probe_chan_ind = find(strcmp(mfile_csd.chan_names, sample_channames{1}));
abs_chan_ind = mfile_csd.chans(1, probe_chan_ind) + chan_offset;
chan_struct = struct(sample_probe, abs_chan_ind);

%% Raw LFP

% Load raw data
data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s_%s.mat', sample_day, sample_time)));
srate = data_mfile.finalSampR;
data_range = timerange_offset(1)*srate:timerange_offset(2)*srate;
lfp_segment = organize_lfp(data_mfile, chan_struct, [], [], [], data_range);
xaxis_lfp = linspace(0, sample_timerange(2) - sample_timerange(1), length(lfp_segment));

h_lfp = figure('Position', figpos);
h_lfp.Renderer = 'painters';
plot(xaxis_lfp, 1000 * lfp_segment, 'k', 'LineWidth', 1);
title('Raw LFP');
xlabel('Time (s)');
ylabel('Potential (uV)');
set(gca, 'FontSize', 16, 'FontName', 'Arial', 'XGrid', 'on', 'Box', 'off');

%% Raw LFP callouts (pre- and post-transition)

mid_sample = 1 + (data_range(end)-data_range(1)) / 2;
callout_offset = srate * 10; % start 10 seconds before and after transition
callout_length = srate * 6; % 6 seconds (1 window)
pre_samps = mid_sample - (callout_offset + (callout_length:-1:1));
post_samps = mid_sample + callout_offset + (1:callout_length);

% plot patches in LFP figure corresponding to callouts
figure(h_lfp);
lfp_ylim = get(gca, 'YLim');
hold on;
x_pre = xaxis_lfp(pre_samps([1, end, end, 1]));
patch(x_pre, repelem(lfp_ylim, 2), repelem(-1, 4), 'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
x_post = xaxis_lfp(post_samps([1, end, end, 1]));
patch(x_post, repelem(lfp_ylim, 2), repelem(-1, 4), 'b', 'FaceAlpha', 0.3, 'LineStyle', 'none');

callout_pos = [0, 0, 295, 84];
h_pre = figure('Position', callout_pos);
plot(lfp_segment(pre_samps), 'k', 'LineWidth', 1);
xticks([]);
yticks([]);
set(gca, 'XColor', 'r');
set(gca, 'YColor', 'r');

h_post = figure('Position', callout_pos);
plot(lfp_segment(post_samps), 'k', 'LineWidth', 1);
xticks([]);
yticks([]);
set(gca, 'XColor', 'b');
set(gca, 'YColor', 'b');

%% Raw spectrogram (unused)

h_raw = figure('Position', figpos);
plot_multitaper(mt_res, struct(...
    'pxx_name',     'pxx', ...
    'take_log',     true, ...
    'xlim',         timerange_offset, ...
    'chan_names',   {sample_fullchannames(1)}, ...
    'axes',         gca));

set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Position', [0.1,0.05,0.75,0.795]);
title('Raw spectrogram');
cb = colorbar;
cb.Label.String = 'Power (mV^2/Hz, dB)';
cb.Position = [0.86, 0.05, 0.03, 0.795];

%% Normalized spectrogram

h_normalized = figure('Position', figpos);
plot_multitaper(mt_res, struct(...
    'pxx_name',     'pxx_rankord', ...
    'xlim',         timerange_offset, ...
    'clim',         [0, 1], ...
    'chan_names',   {sample_fullchannames(1)}, ...
    'axes',         gca));

set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Position', [0.1,0.05,0.75,0.795]);
title('Normalized spectrogram');
cb = colorbar;
cb.Label.String = 'Power rank';
cb.Position = [0.86, 0.05, 0.03, 0.795];

%% NMF loadings and scores

nmf_mfile = matfile(fullfile(sr_dirs.results, sample_day, 'nmf_res.mat'));
sample_chaninds = cellfun(@(cn) find(strcmp(cn, nmf_mfile.chan_names)), sample_fullchannames);

loadings = nmf_mfile.nmf_U;
loadings = loadings{1}{sample_chaninds(1)};
h_loadings = figure('Position', [0, 0, 218, 238]);
sanePColor(1:size(loadings, 2), nmf_mfile.freq_axis, loadings, false, true);
set(gca, 'YScale', 'log', 'FontSize', 16, 'FontName', 'Arial');
xticks(1:size(loadings, 2));
xlabel('Component #');
ylabel('Frequency (Hz)');
title('NMF loadings');
yticks([1, 10, 100]);
box off;

% Loadings for both channels (Figure 4A and first channel is Figure 3D)
all_hr_chan_names = util.make_hr_chan_names(nmf_mfile.chan_names, mt_res.chan_locs);
depths = all_hr_chan_names(sample_chaninds);
scores = nmf_mfile.nmf_V;
scores = scores{1}(sample_chaninds);
n_chans = length(sample_chaninds);
h_scores = gobjects(n_chans, 1);
b_time = nmf_mfile.time_axis >= sample_timerange(1) & nmf_mfile.time_axis <= sample_timerange(2);
score_submats = cellfun(@(sc) sc(b_time, :).', scores, 'uni', false);

for kC = 1:n_chans
    h_scores(kC) = figure('Position', figpos);
    n_comps = size(score_submats{kC}, 1);
    sanePColor(xaxis, 1:n_comps, score_submats{kC});
    title(sprintf('%s, depth = %s um', sample_region, depths(kC)));
    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Position', [0.1,0.05,0.75,0.795]);
    yticks(1:n_comps);
    ylabel('Component #');
    xlabel('Time (s)');
    cb = colorbar;
    cb.Label.String = 'Component score';
    cb.Position = [0.86, 0.05, 0.03, 0.795];
    box off;
end

%% Reconstruction

recon = loadings * score_submats{1};
h_recon = figure('Position', figpos);
sanePColor(xaxis, nmf_mfile.freq_axis, recon, false, true);
set(gca, 'YScale', 'log', 'FontSize', 16, 'FontName', 'Arial', 'Position', [0.1, 0.05, 0.75, 0.795]);
title('NMF reconstruction');
yticks([1, 10, 100]);
ylabel('Frequency (Hz)');
xlabel('Time (s)');
cb = colorbar;
cb.Label.String = 'Power rank';
cb.Position = [0.86, 0.05, 0.03, 0.795];
box off;

%% Discrete state plot

classes = nmf_mfile.filtered_classes;
classes = classes{1};
classes = cellfun(@(cls) cls(b_time), classes, 'uni', false);

h_state = figure('Position', figpos);
class_plot(xaxis, classes);
legend('off');
set(gca, 'FontSize', 16, 'FontName', 'Arial', 'Position', [0.1, 0.05, 0.75, 0.795]);
title('Discrete states');
ylabel('Channel');

%% Discrete state transitions

% open the figure saved during full_analysis_byday since the transition table isn't saved anywhere
h_trans = openfig(fullfile(sr_dirs.results, sample_day, ['transitions_w_sync_', sample_day, '.fig']));
h_trans.Position = figpos;
ax = gca;
set(ax, 'FontSize', 16, 'FontName', 'Arial', 'Position', [0.1, 0.05, 0.75, 0.795]);
cb = h_trans.Children(1);
cb.Position = [0.86, 0.05, 0.03, 0.795];

% Make figure bigger without changing the size of its contents (to fit x axis)
ax.Units = 'pixels';
cb.Units = 'pixels';
shift_px = 50;
h_trans.Position(4) = h_trans.Position(4) + shift_px;
ax.Position(2) = ax.Position(2) + shift_px;
cb.Position(2) = cb.Position(2) + shift_px;
ax.Units = 'relative';
cb.Units = 'relative';

% Zoom in on the relevant segment
xlim(sample_timerange);
xticks(linspace(sample_timerange(1), sample_timerange(2), 7));
xticklabels(linspace(0, sample_timerange(2) - sample_timerange(1), 7));

% Y axis
n_chans = length(classes);
ylim([0, n_chans + 1]);
yticks(1:n_chans);
yticklabels([]);
grid on;

% make lines thicker
set(ax.Children, 'Linewidth', 1.5);

title('Discrete state transitions');
ylabel('Channel');
