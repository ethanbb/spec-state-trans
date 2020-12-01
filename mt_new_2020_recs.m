% Do CSD and multitaper on the October 2020 recordings

sr_dirs = prepSR;

rec_names = {
    '2020-10-26_11-30-00'
    '2020-10-26_13-31-00'
    '2020-10-26_15-33-00'
    '2020-10-27_13-04-00'
    '2020-10-27_15-05-00'
%     '2020-10-28_12-29-00'
    '2020-10-28_14-31-00'
    '2020-10-28_16-45-00'
    '2020-10-29_12-05-00'
    '2020-10-29_14-06-00'
%     '2020-10-29_16-07-00'
    };
n_recs = length(rec_names);

temp = split(rec_names, '_');
recs_dates_times = [rec_names, temp];
rec_dates = unique(temp(:, 1));

%% Do CSD

probe_s = struct('Probe1', 'M1', 'Probe2', 'V1');
for kD = 1:length(rec_dates)
    plot_csd(rec_dates{kD}, probe_s);
end

%% Do multitaper (based on results of pick_channels_from_csd.m)

rec_mt_info = table(repmat({zeros(0, 2)}, n_recs, 1), cell(n_recs, 1), ...
    struct('Probe1', cell(n_recs, 1), 'Probe2', cell(n_recs, 1)), ...
    recs_dates_times(:, 2), recs_dates_times(:, 3), ...
    'VariableNames', {'artifacts', 'chan_names', 'chans', 'date', 'time'}, ...
    'RowNames', rec_names);

% channel notes for each day:
% 10/26 ch. 93 (V1 29) looks somewhat flat/broken
% 10/27 some HF stuff on ch. 101-106 (V1 37-42), could be noise
% 10/28 similar to 10/26 for ch. 93 (V1 29)
% 10/29 looks good.
% fill in channels based on CSD (at least initially)
for kD = 1:length(rec_dates)
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
    if strcmp(date, '2020-10-27')
        % Use a different channel for M1 L5 - the one picked by CSD has some odd artifacts.
        for rec_i = find(b_rec(:)')
            rec_mt_info.chans(rec_i).Probe1(strcmp(csd_info_M1.chan_names, 'L5')) = 33;
        end
    end
end

rec_mt_info.artifacts{'2020-10-27_13-04-00'} = [
    16.8, 20
    135, 156 % "black swan"
    ];

rec_mt_info.artifacts{'2020-10-29_14-06-00'} = [
    1879, 1883
    3705, 3757
    ];

%% Loop through recordings
for kR = 1:length(rec_names)
    %% Do low-resolution analysis first
    rec_name = rec_names{kR};
    data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec_name)));
    save_dir = fullfile(sr_dirs.results, rec_mt_info.date{rec_name}, rec_mt_info.time{rec_name});
    
    options = struct;
    options.artifacts = rec_mt_info.artifacts{rec_name};
    options.chan_names = rec_mt_info.chan_names{rec_name};
    options.chans = rec_mt_info.chans(rec_name);
    options.save = false;
    
    mt_res_lores = multitaper_analysis(data_mfile, options);
    
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
    options.savedir = save_dir;
    options.filename = 'mt_res_layers.mat';
    mt_res = multitaper_analysis(data_mfile, options);
    
end
