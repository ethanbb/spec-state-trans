% Do CSD and multitaper on the 2021 bilateral V1 recordings -
% both baseline and stimulated.

sr_dirs = prepSR;

base_rec_names = {
%     '2021-01-25_13-30-00'  % exclude 1/25 - left side is way too shallow or not in V1
%     '2021-01-27_15-05-00'  % not yet preprocessed
    '2021-01-29_13-46-00'
    '2021-01-29_15-22-00'
    '2021-01-31_14-12-00'
    '2021-01-31_15-41-00'
    '2021-02-01_14-49-00'
    };
n_base_recs = length(base_rec_names);

stim_rec_names = {
%     '2021-01-25_14-32-00'
%     '2021-01-25_16-02-00'
    '2021-01-27_11-00-00'
    '2021-01-27_12-32-00'
    '2021-01-27_14-02-00'  % no snippits currently
    '2021-01-29_10-39-00'
    '2021-01-29_12-09-00'
    '2021-01-29_16-53-00'
    '2021-01-31_09-52-00'
    '2021-01-31_11-21-00'
    '2021-01-31_12-39-00'
    '2021-02-01_12-10-00'
    '2021-02-01_13-48-00'
    '2021-02-02_10-42-00'
    '2021-02-02_12-13-00'
    '2021-02-02_13-43-00'
    '2021-02-02_15-14-00'
    };
n_stim_recs = length(stim_rec_names);

rec_names = sort([base_rec_names; stim_rec_names]);
n_recs = n_base_recs + n_stim_recs;

recs_dates_times = [rec_names, split(rec_names, '_')];
rec_dates = unique(recs_dates_times(:, 2));

%% Do CSD

probe_s = struct('Probe1', 'V1L', 'Probe2', 'V1R');

dead_chans_s = [
%     struct('Probe1', [39, 40, 53], 'Probe2', [23, 24, 30, 31]) % 2021-01-25
    struct('Probe1', [], 'Probe2', [39, 53])      % 2021-01-27
    struct('Probe1', [], 'Probe2', [39, 40, 53])  % 2020-01-29
    struct('Probe1', [], 'Probe2', [39, 53])      % 2020-01-31
    struct('Probe1', [], 'Probe2', [39, 53])      % 2020-02-01
    struct('Probe1', 29, 'Probe2', [])            % 2020-02-02
    ];

for kD = 1:length(rec_dates)
    plot_csd(rec_dates{kD}, probe_s, dead_chans_s(kD));
end

%% Pick channels - L4 and steps of 140 um up and down

depths_um = 140 * (-8:8);
layer_names = [
    arrayfun(@(k) ['Sup', num2str(k)], 8:-1:1, 'uni', false), {'L4'}, ...
    arrayfun(@(k) ['Inf', num2str(k)], 1:8, 'uni', false)
    ];

for kD = 1:length(rec_dates)
    rec_dir = fullfile(sr_dirs.results, rec_dates{kD});
    for side = 'LR'
        pick_csd_channels(rec_dir, depths_um, layer_names, ['V1', side], {}, true);
    end
end

%% Scan for burst suppression

bs_artifacts = cell(n_recs, 1);

%%
for kR = 1:n_recs
    %%
    data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec_names{kR})));
    n = size(data_mfile, 'meanSubFullTrace', 1);
    bs_dur_thresh = 2;
    bs_segments = find_likely_bs(data_mfile, round(linspace(1, n, 4)), bs_dur_thresh);
    bs_artifacts{kR} = bs_segments;

    if exist('check_bs', 'var') && check_bs
        disp(bs_segments);
        % plot using eeglab to check
        lfp = organize_lfp(data_mfile, 1:8:n);
        eegplot(lfp, 'srate', 1000, 'winlength', 20);
    end
end
