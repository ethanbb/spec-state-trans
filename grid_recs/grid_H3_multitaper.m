% Do multitaper analysis on each long probe/grid dataset, using every other channel.
% Start with lo-res version in order to identify any noise channels or artifacts.

%% Prepare
sr_dirs = prepSR;

% common parameters
save_dir = fullfile(sr_dirs.results, 'grid_recs', 'mt_res');

fork_chans = 1:7:64;

grid_chans_mod = reshape(mod(1:66, 7), 11, 6); % includes 2 "dummy" channels
b_dummy_grid_chans = [false(10, 6); false, false, true, true, false, false];
grid_chans_mod = grid_chans_mod(~b_dummy_grid_chans);
grid_chans = find(grid_chans_mod == 1);

default_chans = struct('G', grid_chans, 'P', fork_chans);

% Make struct array of recs to process
% These recs are the baseline ones that don't have BS

rec_names = {
    % Not using 03-07 b/c too much burst suppresion
    '2019-04-25_13-15-00'
    '2019-04-25_14-29-00'
    '2019-08-27_14-52-00'
    };
n_recs = length(rec_names);

grid_recs = table(repmat({zeros(0, 2)}, n_recs, 1), repmat(default_chans, n_recs, 1), ... 
    'VariableNames', {'artifacts', 'chans'}, 'RowNames', rec_names);

rec = '2019-04-25_13-15-00';
grid_recs.artifacts{rec} = [
    454, 1121
    ];
% to avoid artifact channels
grid_recs.chans(rec).G = find([
    grid_chans_mod(1:32) == 4
    grid_chans_mod(33:end) == 0
    ]);

rec = '2019-04-25_14-29-00';
grid_recs.artifacts{rec} = [
    109, 399
    567, -1 % to end
    ];
grid_recs.chans(rec).G = find([
    grid_chans_mod(1:32) == 4
    grid_chans_mod(33:end) == 0
    ]);

rec = '2019-08-27_14-52-00';
grid_recs.artifacts{rec} = [
    0, 511
    1095, 1464
    1806, 1813
    1933, 2316
    2996, 3425
    3749, 3750
    4123, 4363
    5211, 5213
    5722, 5837
    6084, 6086
    6711, 6713
    6984, 7194
    7698, -1
    ];
% designated "noise channels": 7, 8, 26, 33, 44, 54, 55, 56, 57
% also noisy in this dataset: 9, 10, 43, 45, 61, 62, 63
grid_recs.chans(rec).G = [4, 11, 18, 25, 32, 38, 46, 51, 59, 64];

%% Loop through recordings
for kR = 1:length(grid_recs)
    %% Prepare dataset
    rec_name = rec_names{kR};
    data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec_name)));
    data_info = data_mfile.info;
    
    %% Do low-resolution analysis first
    
    options = struct;

    options.artifacts = grid_recs.artifacts{rec_name};

    % get struct of all channels, not including noise chans in the info struct, to figure out how to
    % proceed
    [~, good_chans] = organize_lfp(data_mfile, 1:data_info.channels, data_info.noiseChannels, false);
    
    % make chan_names based on the actual channels going to be used
    % figure out which one is the grid
    probe_names = fieldnames(good_chans);
    options.chan_names = cell(2, 1);
    options.chans = struct;
    for kP = 1:2
        pname = probe_names{kP};
        if startsWith(data_info.([pname, 'Name']), 'E64')
            name_prefix = 'G';
        else
            name_prefix = 'P';
        end
        
        if ~all(ismember(grid_recs.chans(rec_name).(name_prefix), good_chans.(pname)))
            keyboard; % figure out how to choose different channels to avoid bad ones
        end
        options.chans.(pname) = grid_recs.chans(rec_name).(name_prefix)(:);
        options.chan_names{kP} = arrayfun(@(c) sprintf('%s%d', name_prefix, c), ...
            options.chans.(pname), 'uni', false);
    end
    options.chan_names = vertcat(options.chan_names{:});
    
    options.save = false;
    
    mt_res_lores = multitaper_analysis(data_mfile, options);
    
%     % (example code to inspect results - modify as necessary)
%     plot_options = struct;
%     plot_options.pxx_name = 'pxx';
%     plot_options.take_log = true;
%     plot_options.chans = 1:6;
%     
%     plot_multitaper(mt_res_lores, plot_options);
    
    %% Do high-res analysis & save
    
    options.window = 6;
    options.padbase = 60;
    options.winstep = 0.1;
    
    options.save = true;
    options.savedir = save_dir;
    options.filename = sprintf('mt_res_%s.mat', rec_name);
    
    mt_res = multitaper_analysis(data_mfile, options);
end
