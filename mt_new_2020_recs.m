% Do CSD and multitaper on the October 2020 recordings

sr_dirs = prepSR;

rec_names = {
    '2020-10-26_11-30-00'
    '2020-10-26_13-31-00'
    '2020-10-26_15-33-00'
%     '2020-10-27_13-04-00'
%     '2020-10-27_15-05-00' # exclude due to hard-to-interpret csd
    '2020-10-28_12-29-00'
    '2020-10-28_14-31-00'
    '2020-10-28_16-45-00'
    '2020-10-29_12-05-00'
%     '2020-10-29_14-06-00' # excluded b/c it's almost all artifact
    '2020-10-29_16-07-00'
    };
n_recs = length(rec_names);

temp = split(rec_names, '_');
recs_dates_times = [rec_names, temp];
rec_dates = unique(temp(:, 1));
n_dates = length(rec_dates);

%% Do CSD of LED trials

probe_s = struct('Probe1', 'M1', 'Probe2', 'V1');

% exclude bad channels noted in recording info (based on depth, artifacts)
bad_chans = cell(n_dates, 1);
for kD = 1:n_dates
    date_recs = recs_dates_times(strcmp(recs_dates_times(:, 2), rec_dates{kD}), 1);
    rec1 = date_recs{1};
    rec1_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec1)));
    rec1_info = rec1_mfile.info;
    noise_chans = rec1_info.noiseChannels;
    
    bad_chans{kD}.Probe1 = noise_chans(noise_chans <= 64);
    bad_chans{kD}.Probe2 = noise_chans(noise_chans > 64) - 64;
    plot_csd(rec_dates{kD}, probe_s, bad_chans{kD});
end

%% Pick channels - L4 and steps of 140 um up and down

depths_um = 140 * (-8:8);
layer_names = [
    arrayfun(@(k) ['Sup', num2str(k)], 8:-1:1, 'uni', false), {'L4'}, ...
    arrayfun(@(k) ['Inf', num2str(k)], 1:8, 'uni', false)
    ];

for kD = 1:n_dates
    pick_csd_channels(fullfile(sr_dirs.results, rec_dates{kD}), ...
        depths_um, layer_names, 'V1', {'M1'}, true, bad_chans{kD}.Probe2, {bad_chans{kD}.Probe1});
end

%% Scan for burst suppression

bs_artifacts = cell(n_recs, 1);

%%
for kR = 1:n_recs
    %%
    data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec_names{kR})));
    data_info = data_mfile.info;
    noise_chans = data_info.noiseChannels;
    n = size(data_mfile, 'meanSubFullTrace', 1);
    non_noisechans = setdiff(1:n, noise_chans);
    chans_to_use = non_noisechans(round(linspace(1, length(non_noisechans), 4)));
    bs_dur_thresh = 2;
    bs_segments = find_likely_bs(data_mfile, chans_to_use, bs_dur_thresh);
    %disp(bs_segments)
    bs_artifacts{kR} = bs_segments;

    % plot using eeglab to check
%     lfp = organize_lfp(data_mfile, 1:8:n);
%     eegplot(lfp, 'srate', 1000, 'winlength', 20);
end

%% Make info table and manual adjustments to artifacts & channels

rec_mt_info = table(bs_artifacts, cell(n_recs, 1), ...
    struct('Probe1', cell(n_recs, 1), 'Probe2', cell(n_recs, 1)), ...
    recs_dates_times(:, 2), recs_dates_times(:, 3), ...
    'VariableNames', {'artifacts', 'chan_names', 'chans', 'date', 'time'}, ...
    'RowNames', rec_names);

% channel notes for each day:
% 10/26 ch. 93 (V1 29) looks somewhat flat/broken
% 10/27 chs. 15 and 93 look bad. some HF stuff on ch. 101-106 (V1 37-42), could be noise
% 10/28 similar to 10/26 for ch. 93 (V1 29). On 3rd rec, 80, 97, 114 (V1 16, 33, 50) have some weird
% artifacts.
% 10/29 looks good.
% fill in channels based on CSD (at least initially)
for kD = 1:n_dates
    date = rec_dates{kD};
    b_rec = strcmp(rec_mt_info.date, date);
    
    % for these recordings, we know Probe1 = M1 and Probe2 = V1
    csd_dir = fullfile(sr_dirs.results, date);
    csd_info_M1 = matfile(fullfile(csd_dir, 'csd_M1.mat'));
    csd_info_V1 = matfile(fullfile(csd_dir, 'csd_V1.mat'));
    
    rec_mt_info.chan_names(b_rec) = {[
        strcat('M1_', csd_info_M1.chan_names), ...
        strcat('V1_', csd_info_V1.chan_names)
        ]};
    rec_mt_info.chans(b_rec) = struct('Probe1', csd_info_M1.chans, 'Probe2', csd_info_V1.chans);
    if strcmp(date, '2020-10-26')
        % Avoid M1 channel 32
        for rec_i = find(b_rec(:)')
            rec_mt_info.chans(rec_i).Probe1(csd_info_M1.chans == 32) = 33;
        end
    end
    
    if strcmp(date, '2020-10-27')
        % Skip V1 Inf4 due to artifacts
        for rec_i = find(b_rec(:)')
            rec_mt_info.chans(rec_i).Probe2(strcmp(csd_info_V1.chan_names, 'Inf4')) = [];
            rec_mt_info.chan_names{rec_i}(strcmp(rec_mt_info.chan_names{rec_i}, 'V1_Inf4')) = [];
        end
    end
    
    if strcmp(date, '2020-10-28')
        % Avoid V1 channel 33
        for rec_i = find(b_rec(:)')
            rec_mt_info.chans(rec_i).Probe2(csd_info_V1.chans == 33) = 32;
        end
    end
end

% artifact changes/additions (aside from detected BS)

% '2020-10-27_13-04-00' - [135, 156] is "black swan", but it's surrounded by burst suppression

rec_mt_info.artifacts{'2020-10-27_13-04-00'} = [
    rec_mt_info.artifacts{'2020-10-27_13-04-00'}
    2326, 2327 % for channel 97
    ];

rec_mt_info.artifacts{'2020-10-27_15-05-00'} = [
    rec_mt_info.artifacts{'2020-10-27_15-05-00'}
    544.5, 546 % for channel 40
    ];

rec_mt_info.artifacts{'2020-10-28_16-45-00'} = [
    rec_mt_info.artifacts{'2020-10-28_16-45-00'}
    703, 705 % for channel 104
    ];
    
% rec_mt_info.artifacts{'2020-10-29_14-06-00'} = [
%     rec_mt_info.artifacts{'2020-10-29_14-06-00'}
%     1879, 1883
%     3705, 3757
%     ];

rec_mt_info.artifacts{'2020-10-29_16-07-00'} = [
    rec_mt_info.artifacts{'2020-10-29_16-07-00'}
    1353, 1360  % cross-channel
    2862, 2867  % cross-channel
    5961, 5963  % cross-channel
    6038, 6040  % for channel 33
    7004, 7017  % cross-channel
    ];


%% Finally do multitaper - loop through recordings
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
