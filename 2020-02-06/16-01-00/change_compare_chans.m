% Do PCA separately for V1 and motor, then get the change speed in each
% using a small # of components and quantify order of change peaks.

%% Setup

redo_all = false; % change to true to recompute from scratch

redo_pca = redo_all || false; % change to true to redo PCA
redo_comps = redo_all || false; % change to true to re-visualize PCA and select components
redo_change = redo_all || false; % change to true to recalculate change speeds
redo_extrema = redo_all || false; % change to true to redo extrema finding

prepSR;

recdate = '2020-02-06';
time = '16-01-00';
res_file = fullfile(results_dir, recdate, time, 'mt_res.mat');
res_mfile = matfile(res_file, 'Writable', true);

mt_opts = res_mfile.options; % necessary due to MatFile quirk
Fw = 1 / mt_opts.winstep; % "window rate"

chans = res_mfile.name;
chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false);
n_chans = numel(chans);

%% Do PCA, just taking 10 components b/c going to decide which ones to keep visually

pc_data = cell(n_chans, 1);

total_comps = 10;

for kC = 1:n_chans
    res_name = sprintf('pxx_pca_1rec_%s', chan_vnames{kC});
    
    if ~redo_pca && isprop(res_mfile, res_name) % just load if possible
        pc_data(kC) = res_mfile.(res_name)(kC, 1);
    else
        pca_opts = struct;
        pca_opts.name = res_name;
        pca_opts.chans = kC;
        pca_opts.thresh_type = 'comps';
        pca_opts.thresh = total_comps;

        [mt_pca_output, fh] = mt_pca(res_file, pca_opts);
        pc_data(kC) = mt_pca_output(kC);

        savefig(fh, sprintf('pca_%s.fig', chan_vnames{kC}));
    end
end

%% Visualize to identify which components to keep (start with scatter 1-3, can do more if necessary)

comp_fname = 'change_compare_comps';

if ~redo_comps && isprop(res_mfile, comp_fname)
    comps2use = res_mfile.(comp_fname);
else
    % these times of interest have been painstakingly manually identified
    hstarts = [1521, 2197, 3722, 4530, 5180, 6081, 7619]; % in seconds
    hends = [1614, 2634, 3814, 4702, 5259, 6240, 7712];   % ditto
    
    Fs_mt = 10; % since our window step is 0.1 seconds, sample freq = 10 Hz
    htimes = arrayfun(@(s, e) Fs_mt*s+1:Fs_mt*e, hstarts, hends, 'uni', false);
    htimes = cell2mat(htimes);
    
    
    comps2use = false(n_chans, total_comps);
    
    figure;
    
    for kC = 1:n_chans
        comps2show = 1:3;
        
        while ~isempty(comps2show)
            
            hold off;
            dscatter(pc_data{kC}(comps2show, :));
            set(gca, 'Interactions', [zoomInteraction, rotateInteraction, rulerPanInteraction]);
            hold on;
            xlabel(sprintf('PC%d', comps2show(1)));
            ylabel(sprintf('PC%d', comps2show(2)));
            
            % show times of interest
            hdata = arrayfun(@(c) pc_data{kC}(c, htimes), comps2show, 'uni', false);
            if length(comps2show) == 3
                zlabel(sprintf('PC%d', comps2show(3)));
                scatter3(hdata{:}, 20, 'r', 'o');
            else
                scatter(hdata{:}, 20, 'r', 'o');
            end
            
            title(sprintf('Principal components in region: %s', chans{kC}));
            
            % Ask user which components to use in analysis
            comp_prompt = 'Enter array of most informative components: ';
            comp_validator = @(comps) all(ismember(comps, 1:total_comps));
            comp_errmsg = sprintf('Must be an array of component numbers from 1 to %d', total_comps);
            comps2add = safeinput(comp_prompt, comp_validator, comp_errmsg);

            comps2use(kC, comps2add) = true;
            
            % Ask user which components to show next
            next_prompt = 'Enter 2 or 3 components to show next, or nothing to move on: ';
            next_validator = @(comps) comp_validator(comps) && ismember(numel(comps), [0, 2, 3]);
            next_errmsg = sprintf('Must be empty or an array of 2-3 component numbers from 1 to %d', total_comps);
            comps2show = safeinput(next_prompt, next_validator, next_errmsg);
        end        
    end
    
    % convert to n_regions x 1 cell
    comps2use = arrayfun(@(kC) find(comps2use(kC, :)), (1:n_chans)', 'uni', false);
    res_mfile.(comp_fname) = comps2use;
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
    change_opts.smooth_method = 'gaussian';
    change_opts.name = change_fname;
    change_opts.savefigs = false;
    
    % more smoothing:
    change_opts.smooth_span = 100;

    % exponential smoothing:
    change_opts.smooth_method = 'exp';
    change_opts.diff_step = 1/Fw; % detect sharp change right at sample of interest

    change_opts.norm_type = 'none'; % get velocity instead of speed (meaningful peaks & troughs)
    
    [change_vel{kC}, change_time] = pca_change(res_file, change_opts);
    title(sprintf('%s PC1 velocity', chans{kC}));
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
    
    ispeak = cellfun(@(p) change_time == p, peaks, 'uni', false);
    istrough = cellfun(@(t) change_time == t, troughs, 'uni', false);
else
    
    ispeak = cell(n_chans, 1);
    istrough = cell(n_chans, 1);
    
    for kC = 1:n_chans
        data_range = max(change_vel{kC}) - min(change_vel{kC});
        ispeak{kC} = islocalmax(change_vel{kC}, 'MinProminence', 0.4 * data_range);
        istrough{kC} = islocalmin(change_vel{kC}, 'MinProminence', 0.4 * data_range);
        
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
