% Do PCA separately for V1 and motor, then get the change speed in each
% using a small # of components and quantify order of change peaks.

%% Setup

redo_all = true; % change to true to recompute from scratch

redo_smooth = redo_all || false; % change to true to redo pre-PCA smoothing
redo_pca = redo_all || false; % change to true to redo PCA
redo_comps = redo_all || false; % change to true to re-visualize PCA and select components
redo_change = redo_all || false; % change to true to recalculate change speeds
redo_extrema = redo_all || false; % change to true to redo extrema finding
redo_match = redo_all || false; % change to true to redo extrema matching

prepSR;

recdate = '2020-03-11';
time = '12-31-00';
res_file = fullfile(results_dir, recdate, time, 'mt_res.mat');
res_mfile = matfile(res_file, 'Writable', true);

mt_opts = res_mfile.options; % necessary due to MatFile quirk
Fw = 1 / mt_opts.winstep; % "window rate"

chans = [2, 5];
chan_names = {'V1', 'MC'};
chan_vnames = cellfun(@matlab.lang.makeValidName, chan_names, 'uni', false);
n_chans = numel(chans);

%% Smooth pxx before doing PCA

smooth_fname = 'pxx_smooth';

if redo_smooth || ~isprop(res_mfile, smooth_fname)
    
    pxx = res_mfile.pxx;
    pxx_smooth = cell(size(pxx));

    for kC = 1:n_chans
        chan_pxx = pxx{chans(kC)};

        % smooth across frequencies with a median filter
        pxx_freqsmooth = medfilt1(chan_pxx, 40);

        % exponential smoothing in time:
        smooth_span = 60; % seconds
        sm_span_samp = smooth_span * Fw;
        decay_2t = normpdf(2.5) / normpdf(0);
        [b_exp, a_exp] = exp_filter(sm_span_samp, decay_2t);

        % manually filter forward and backward since filtfilt can't deal with nans or matrices
        pxx_smooth{chans(kC)} = filtfilt_segs(b_exp, a_exp, pxx_freqsmooth, res_mfile.seg_windows, 2);

        % gaussian smoothing
%         pxx_smooth{chans(kC)} = smoothdata(pxx_freqsmooth, 2, 'gaussian', sm_span_samp, 'includenan');
    end

    res_mfile.(smooth_fname) = pxx_smooth;
end

%% Plot and save smoothed data

plot_opts = struct;
plot_opts.pxx_name = smooth_fname;
plot_opts.filename = 'multitaper_sm.fig';
plot_opts.chans = chans;

plot_multitaper(res_file, plot_opts);

%% Concatenate related channels along frequency axis

% res_mfile.pxx_smooth_cat = {
%     cell2mat(res_mfile.(smooth_fname)(1:3, 1))
%     cell2mat(res_mfile.(smooth_fname)(4:6, 1))
%     };
%
% chan_names = {'V1'; 'MC'};
% chan_vnames = chan_names;
% n_chans = numel(chan_names);

%% Do PCA, just taking 10 components b/c going to decide which ones to keep visually

pc_data = cell(n_chans, 1);

total_comps = 10;

for kC = 1:n_chans
    res_name = sprintf('pxx_pca_1rec_%s', chan_vnames{kC});
    
    if ~redo_pca && isprop(res_mfile, res_name) % just load if possible
        pc_data(kC) = res_mfile.(res_name)(chans(kC), 1);
    else
        pca_opts = struct;
        pca_opts.pxx_name = smooth_fname;
        pca_opts.name = res_name;
        pca_opts.chans = chans(kC);
        pca_opts.thresh_type = 'comps';
        pca_opts.thresh = total_comps;

        [pc_data(kC), fh] = mt_pca(res_file, pca_opts);

        savefig(fh, sprintf('pca_%s.fig', chan_vnames{kC}));
    end
end

%% Actually just use the component with highest variance:
comps2use = repmat({1}, n_chans, 1);

%% Compute change velocity for each channel using pca_change

change_vel = cell(n_chans, 1);
change_fhs = gobjects(n_chans, 1); % figure handles

for kC = 1:n_chans
    change_fname = sprintf('pca_change_%s', chan_vnames{kC});
    change_figname = fullfile(results_dir, recdate, time, [change_fname, '.fig']);
    
    if ~redo_change && isprop(res_mfile, change_fname) && exist(change_figname, 'file')
        change_vel{kC} = res_mfile.(change_fname);
        change_time = res_mfile.([change_fname, '_time']);
        change_fhs(kC) = openfig(change_figname);
        continue;
    end
    
    change_opts = struct;
    change_opts.pca_name = sprintf('pxx_pca_1rec_%s', chan_vnames{kC});
    change_opts.comps = comps2use{kC};
    change_opts.name = change_fname;
    change_opts.savefigs = false;
    
    % no smoothing here:
    change_opts.smooth_span = 1;

    % regular diff
    change_opts.diff_step = 1/Fw; % detect sharp change right at sample of interest

    change_opts.norm_type = 'none'; % get velocity instead of speed (meaningful peaks & troughs)
    
    [change_vel{kC}, change_time] = pca_change(res_file, change_opts);
    title(sprintf('%s PC1 velocity', chan_names{kC}));
    change_fhs(kC) = gcf;
    savefig(change_fhs(kC), change_figname);
end

%% Sanity check - plot over spectrogram

if ~exist('h_mt', 'var') || ~isvalid(h_mt)
    h_mt = openfig(fullfile(results_dir, recdate, time, 'multitaper.fig'));
end

% deal with figures that don't include unnormalized plots
n_cols = numel(h_mt.Children) / n_chans;

for kC = 1:n_chans
    figure(h_mt);
    hax = subplot(n_chans, n_cols, n_cols*kC);

    % clear except for surface plot
    ax_plots = allchild(hax);
    delete(ax_plots(1:end-1));
    hold on;

    yyaxis right;
    plot(change_time, change_vel{kC}, 'k-');
    ylim([min(0, min(change_vel{kC})), 2 * max(change_vel{kC})]);
    ylabel('Change per second');
end

%% Invert change velocity of one or more channels if necessary at this point.
% (depending on sign of first PCA component)
% peaks should correspond to transitions to state with higher relative power in higher freqs
% and vice versa.

if redo_change
    % Manually edit below if necessary:
    chans2invert = [];
    for kC = chans2invert
        change_vel{kC} = -change_vel{kC};

        change_fname = sprintf('pca_change_%s', chan_vnames{kC});
        res_mfile.(change_fname) = -res_mfile.(change_fname);
    end
end

%% Find change extrema

peak_fname = 'vel_peaks';
trough_fname = 'vel_troughs';

% 'peaks' and 'troughs' saved in file are the times of change extrema in seconds.
% 'ispeak' and 'istrough' are logical arrays marking these events, on the time scale of change_vel.

if ~redo_extrema && isprop(res_mfile, peak_fname) && isprop(res_mfile, trough_fname)
    peaks = res_mfile.(peak_fname);
    troughs = res_mfile.(trough_fname);
    
    ispeak = cellfun(@(p) ismember(change_time, p), peaks, 'uni', false);
    istrough = cellfun(@(t) ismember(change_time, t), troughs, 'uni', false);
else
    
    ispeak = cell(n_chans, 1);
    istrough = cell(n_chans, 1);
    
    for kC = 1:n_chans
        vel_iqr = iqr(change_vel{kC});
        ispeak{kC} = islocalmax(change_vel{kC}, 'MinProminence', 3.75 * vel_iqr);
        istrough{kC} = islocalmin(change_vel{kC}, 'MinProminence', 3.75 * vel_iqr);
    end
    
    peaks = cellfun(@(bp) change_time(bp), ispeak, 'uni', false);
    troughs = cellfun(@(bt) change_time(bt), istrough, 'uni', false);
    
    res_mfile.(peak_fname) = peaks;
    res_mfile.(trough_fname) = troughs;
end

% plot extrema for each channel
for kC = 1:n_chans
    figure(change_fhs(kC));
    
    % remove previous extrema if necessary
    ax_plots = allchild(gca);
    delete(ax_plots(1:end-1));    
    
    hold on;
    plot(peaks{kC}, change_vel{kC}(ispeak{kC}), 'ro', 'MarkerSize', 10);
    plot(troughs{kC}, change_vel{kC}(istrough{kC}), 'mo', 'MarkerSize', 10);   
end

%% Match peaks and troughs between channels

% Algorithm:
% * for each peak and trough, find equivalent in other channel that is the closest in time
% * discard any that are more than max_dist away
% * discard any that are not reciprocal
% * the remainders are our matches.

% Here we assume there are 2 channels (TODO if necessary - extend so it works with more)

matched_peaks_fname = 'matched_peaks';
matched_troughs_fname = 'matched_troughs';

max_dist = 30; % dist b/w matching peak and trough in seconds

if ~redo_match && isprop(res_mfile, matched_peaks_fname) && isprop(res_mfile, matched_troughs_fname)
    matched_peaks = res_mfile.(matched_peaks_fname);
    matched_troughs = res_mfile.(matched_troughs_fname);
else

    res_mfile.match_max_dist = max_dist;

    nearest_peak_ind = cell(2, 1);
    nearest_trough_ind = cell(2, 1);
    
    n_peaks = cellfun('length', peaks);
    n_troughs = cellfun('length', troughs);
    
    for kC = 1:2
        nearest_peak_ind{kC} = zeros(1, n_peaks(kC));
        
        for kP = 1:n_peaks(kC)
            [dist, nearest] = min(abs(peaks{kC}(kP) - peaks{3-kC}));
            if ~isempty(nearest) && dist <= max_dist
                nearest_peak_ind{kC}(kP) = nearest;
            end
        end
        
        nearest_trough_ind{kC} = zeros(1, n_troughs(kC));
        
        for kT = 1:n_troughs(kC)
            [dist, nearest] = min(abs(troughs{kC}(kT) - troughs{3-kC}));
            if ~isempty(nearest) && dist <= max_dist
                nearest_trough_ind{kC}(kT) = nearest;
            end
        end
    end
    
    % now that both sides have found their potential matches, loop again
    for kC = 1:2
        % set "nearest peak" entries for which the match is not reciprocal to 0
        % those that are already 0 will see the augmented 0 and thus will not match.
        augmented_matches = [0, nearest_peak_ind{3-kC}];
        nearest_peak_ind{kC}(augmented_matches(nearest_peak_ind{kC} + 1) ~= 1:n_peaks(kC)) = 0;
        
        augmented_matches = [0, nearest_trough_ind{3-kC}];
        nearest_trough_ind{kC}(augmented_matches(nearest_trough_ind{kC} + 1) ~= 1:n_troughs(kC)) = 0;
    end
    
    res_mfile.num_no_match = ...
        sum(cellfun(@(npi) sum(npi==0), nearest_peak_ind)) + ...
        sum(cellfun(@(nti) sum(nti==0), nearest_trough_ind));

    % finally, get the actual times
    n_matched_peaks = sum(nearest_peak_ind{1} > 0);
    n_matched_troughs = sum(nearest_trough_ind{1} > 0);
    
    res_mfile.num_match = n_matched_peaks * 2 + n_matched_troughs * 2;

    matched_peaks = zeros(2, n_matched_peaks);
    matched_troughs = zeros(2, n_matched_troughs);
    
    for kC = 1:2
        matched_peaks(kC, :) = peaks{kC}(nearest_peak_ind{3-kC}(nearest_peak_ind{3-kC} > 0));
        matched_troughs(kC, :) = troughs{kC}(nearest_trough_ind{3-kC}(nearest_trough_ind{3-kC} > 0));
    end
    
    res_mfile.(matched_peaks_fname) = matched_peaks;
    res_mfile.(matched_troughs_fname) = matched_troughs;
end

% Visualize the matched peaks and troughs
figure('Position', [0, 0, 800, 350]);

shifted_vel = cell(2, 1);
shifted_vel{1} = change_vel{1} - min(change_vel{1});
shifted_vel{2} = change_vel{2} - max(change_vel{2});

plot([change_time(1), change_time(end)], [0, 0], 'k');

hold on;

for kC = 1:2
    plot(change_time, shifted_vel{kC}, 'b');    
    plot(matched_peaks(kC, :), shifted_vel{kC}(ismember(change_time, matched_peaks(kC, :))), ...
        'r.', 'MarkerSize', 15);
    plot(matched_troughs(kC, :), shifted_vel{kC}(ismember(change_time, matched_troughs(kC, :))), ...
        'm.', 'MarkerSize', 15);   
end

axis tight;
ylims = get(gca, 'YLim');
ylim([-max(abs(ylims)), max(abs(ylims))]);

% plot connections
for kP = 1:size(matched_peaks, 2)
    plot(matched_peaks(:, kP), ...
        [shifted_vel{1}(change_time == matched_peaks(1, kP)), shifted_vel{2}(change_time == matched_peaks(2, kP))], ...
        'r', 'LineWidth', 2);
end

for kT = 1:size(matched_troughs, 2)
    plot(matched_troughs(:, kT), ...
        [shifted_vel{1}(change_time == matched_troughs(1, kT)), shifted_vel{2}(change_time == matched_troughs(2, kT))], ...
        'm', 'LineWidth', 2);
end

ylabel(sprintf('%s     |     %s', chan_names{2}, chan_names{1}), 'FontSize', 16);

xlabel('Time (s)');
pct_matched = res_mfile.num_match / (res_mfile.num_no_match + res_mfile.num_match) * 100;
title(sprintf('Matched spectral change peaks and troughs, %.0f%% matched (%s_%s)',...
    pct_matched, recdate, time), 'Interpreter', 'none');
savefig('matched_change.fig');