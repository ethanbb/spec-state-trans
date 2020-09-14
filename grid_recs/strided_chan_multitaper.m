% Do multitaper analysis on each long probe/grid dataset, using every other channel.
% Start with lo-res version in order to identify any noise channels or artifacts.

%% Prepare
sr_dirs = prepSR;

% These recs are the baseline ones that don't have BS
grid_recs = {
    '2019-03-07_12-50-00'
%     '2019-03-14_13-52-00'
%     '2019-04-25_13-15-00'
%     '2019-04-25_14-29-00'
%     '2019-08-27_14-52-00'
    };

% common parameters
save_dir = fullfile(sr_dirs.results, 'grid_recs', 'mt_res');
fork_chans = 1:2:64;
b_grid_chans = reshape(mod(1:66, 2) == 1, 11, 6); % includes 2 "dummy" channels
b_dummy_grid_chans = [false(10, 6); false, false, true, true, false, false];
grid_chans = find(b_grid_chans(~b_dummy_grid_chans));

grid_chan_names = arrayfun(@(c) sprintf('G%d', c), grid_chans(:), 'uni', false);
fork_chan_names = arrayfun(@(c) sprintf('F%d', c), fork_chans(:), 'uni', false);

chans = struct('Probe1', grid_chans, 'Probe2', fork_chans);

%% Loop through recordings
for kR = 1:length(grid_recs)
    %% Prepare dataset
    data_mfile = matfile(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', grid_recs{kR})));
    
    %% Do low-resolution analysis first
    
    options = struct;
    options.chans = chans;
    
    options.chan_names = [grid_chan_names; fork_chan_names];
    
    options.save = false;
    
    mt_res_lores = multitaper_analysis(data_mfile, options);
    
    % (example code to inspect results - modify as necessary)
    % plot_options = struct;
    % plot_options.pxx_name = 'pxx';
    % plot_options.take_log = true;
    % plot_options.chans = 1:6;
    % 
    % plot_multitaper(mt_res_lores, plot_options);
    
    %% Do high-res analysis & save
    
    options.window = 6;
    options.padbase = 60;
    options.winstep = 0.1;
    
    options.save = true;
    options.savedir = save_dir;
    options.filename = sprintf('mt_res_%s.mat', grid_recs{kR});
    
    mt_res = multitaper_analysis(data_mfile, options);
end
