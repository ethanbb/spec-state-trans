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
%     '2021-02-01_14-49-00' % exclude 2/1 due to a combination of BS and noisy channels
    };
n_base_recs = length(base_rec_names);

stim_rec_names = {
%     '2021-01-25_14-32-00' % see above
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
%     '2021-02-01_12-10-00' % see above
%     '2021-02-01_13-48-00'
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
n_dates = length(rec_dates);

%% Scan for burst suppression

bs_artifacts = cell(n_recs, 1);

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

rec_mt_info = table(bs_artifacts, cell(n_recs, 1), ...
    struct('Probe1', cell(n_recs, 1), 'Probe2', cell(n_recs, 1)), ...
    recs_dates_times(:, 2), recs_dates_times(:, 3), ...
    'VariableNames', {'artifacts', 'chan_names', 'chans', 'date', 'time'}, ...
    'RowNames', rec_names);

%% Do CSD

probe_s = struct('Probe1', 'V1L', 'Probe2', 'V1R');

bad_chans_t = table(struct('Probe1', cell(n_dates, 1), 'Probe2', cell(n_dates, 1)), ...
    'VariableNames', {'bad_chans'}, 'RowNames', rec_dates);

common_bad_r_chan_dates = {'2021-01-27', '2021-01-29', '2021-01-31'};
for kD = 1:length(common_bad_r_chan_dates)
    bad_chans_t.bad_chans(common_bad_r_chan_dates{kD}).Probe2 = [39, 53];
end

bad_chans_t.bad_chans('2021-01-29').Probe2(end+1) = 40;
bad_chans_t.bad_chans('2021-02-02').Probe1 = 29;

stim_event_chan_s = struct('Probe1', 2, 'Probe2', 1);

for kD = 1:length(rec_dates)
    plot_csd(rec_dates{kD}, probe_s, bad_chans_t.bad_chans(kD), [], rec_mt_info, stim_event_chan_s);
end

%% Pick channels - L4 and steps of 140 um up and down

depths_um = 140 * (-8:8);
layer_names = [
    arrayfun(@(k) ['Sup', num2str(k)], 8:-1:1, 'uni', false), {'L4'}, ...
    arrayfun(@(k) ['Inf', num2str(k)], 1:8, 'uni', false)
    ];

%%
for kD = 1:length(rec_dates)
    %%
    rec_dir = fullfile(sr_dirs.results, rec_dates{kD});
    for side = 'LR'
        pick_csd_channels(rec_dir, depths_um, layer_names, ['V1', side], {}, true);
    end
end

%% Add channel info to multitaper info table

for kD = 1:length(rec_dates)
    rec_date = rec_dates{kD};
    b_rec = strcmp(rec_mt_info.date, rec_date);
    
    % for these recordings, Probe1 = V1L and Probe2 = V1R.
    csd_dir = fullfile(sr_dirs.results, rec_date);
    csd_info_V1L = matfile(fullfile(csd_dir, 'csd_V1L.mat'));
    csd_info_V1R = matfile(fullfile(csd_dir, 'csd_V1R.mat'));
    
    rec_mt_info.chan_names(b_rec) = {[
        strcat('V1L_', csd_info_V1L.chan_names), ...
        strcat('V1R_', csd_info_V1R.chan_names)  ...
        ]};
    rec_mt_info.chans(b_rec) = struct('Probe1', csd_info_V1L.chans, 'Probe2', csd_info_V1R.chans);
    
    % make any additional per-date channel changes here
end

% add any additional artifacts here
rec_mt_info.artifcats{'2021-01-27_11-00-00'} = [
    rec_mt_info.artifacts{'2021-01-27_11-00-00'}
    2248, 2250
    2659, 2661
    3733, 3735
    4107, 4109
    5258, 5259
    ];

rec_mt_info.artifcats{'2021-01-27_12-32-00'} = [
    rec_mt_info.artifacts{'2021-01-27_12-32-00'}
    166,  167
    513,  515
    1327, 1328
    1398, 1399
    3248, 3249
    ];

rec_mt_info.artifacts{'2021-01-31_09-52-00'} = [
    rec_mt_info.artifacts{'2021-01-31_09-52-00'}
    3956, 3958
    ];

rec_mt_info.artifacts{'2021-01-31_15-41-00'} = [
    rec_mt_info.artifacts{'2021-01-31_15-41-00'}
    3950, 3951
    4264, 4266
    ];

%% Loop through recordings and do multitaper
for kR = 1:n_recs
    %% Do low-resolution analysis first
    rec_name = rec_names{kR};
    rec_date = rec_mt_info.date{rec_name};
    rec_time = rec_mt_info.time{rec_name};
    
    data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec_name)));
    
    options = struct;
    options.artifacts = rec_mt_info.artifacts{rec_name};
    options.chan_names = rec_mt_info.chan_names{rec_name};
    options.chans = rec_mt_info.chans(rec_name);
    options.save = false;
 
    %%
%     mt_res_lores = multitaper_analysis(data_mfile, options);
    
%     % (example code to inspect results - modify as necessary)
%     plot_options = struct;
%     plot_options.pxx_name = 'pxx';
%     plot_options.take_log = true;
%     plot_options.chans = 1:4;
%     
%     plot_multitaper(mt_res_lores, plot_options);

    %% Do high-res analysis & save
    
    options.window = 6;
    options.padbase = 60;
    options.winstep = 0.1;
    options.save = true;
    options.savedir = fullfile(sr_dirs.results, rec_date, rec_time);
    options.filename = 'mt_res_layers.mat';
    mt_res = multitaper_analysis(data_mfile, options);

    %% Version with CSDs
    options.window = 6;
    options.padbase = 60;
    options.winstep = 0.1;
    options.save = true;
    options.filename = 'mt_res_layers.mat';
    options.use_csd = true;
    options.bad_chans = bad_chans_t.bad_chans(rec_date);
    options.savedir = fullfile(sr_dirs.results, [rec_date, '_csd'], rec_time);

    raw_mfile = matfile(fullfile(sr_dirs.raw, rec_date, 'matlab', [rec_name, '.mat']));

    multitaper_analysis(raw_mfile, options);
end
