function [mt_file_paths, lfp_inds] = plot_example_data_from_classes(rec_date, channel, plot_length, n_to_plot)

if ~exist('n_to_plot', 'var') || isempty(n_to_plot)
    n_to_plot = 1;
end

sr_dirs = prepSR;

Fs = 1000;

inds_mfile = matfile('characteristic_inds.mat');

[~, input_s_all] = gather_exp_info;
day_ind = find(strcmp(rec_date, inds_mfile.run_dates));
this_input_s = input_s_all(day_ind);
n_files = length(this_input_s.mt_res_in);

% For given channel & day, convert window indices of each clean segment to
% start and ending windows, indexed after noise windows have been removed.
% (This corresponds to sample numbers after multitaper and concatenation.)
segment_info = table(cell(n_files, 1), (1:n_files)', 'VariableNames', {'data', 'filenum'});
total_window_offset = 0;

for kF = 1:n_files    
    mt_mfile = matfile(this_input_s.mt_res_in{kF});
    
    if kF == 1
        opts = mt_mfile.options;
        
        % default to just 1 window
        if ~exist('plot_length', 'var') || isempty(plot_length)
            plot_length = opts.window;
        end
        plot_samples = plot_length * Fs;
        samps_in_window = opts.window * Fs;
        step_samps = opts.winstep * Fs;
        plot_windows = 1 + ceil(max(0, (plot_samples - samps_in_window) / step_samps));
        chan_names = mt_mfile.name;
        chan_ind = find(strcmp(channel, chan_names), 1);
        assert(~isempty(chan_ind), 'Channel %s not found', channel);
        
        % annoying thing to get the struct for requested channel, for later
        cum_chans = cumsum(structfun(@length, opts.chans));
        kProbe = find(chan_ind <= cum_chans, 1);
        pre_chans = [0; cum_chans(:)];
        rel_chan_ind = chan_ind - pre_chans(kProbe);
        probe_names = fieldnames(opts.chans);
        chan_struct = struct(probe_names{kProbe}, opts.chans.(probe_names{kProbe})(rel_chan_ind));
    end
    
    seg_windows = mt_mfile.seg_windows;
    seg_samples = mt_mfile.seg_samples;
    n_segs = length(seg_windows);
    segment_info.data{kF} = table(...
        'Size', [n_segs, 3], ...
        'VariableTypes', repmat({'double'}, 1, 3), ...
        'VariableNames', {'first_window', 'last_window', 'first_window_start_samp'});
                      
    for kS = 1:n_segs
        n_windows = length(seg_windows{kS});
        segment_info.data{kF}.first_window(kS) = total_window_offset + 1;
        segment_info.data{kF}.last_window(kS) = total_window_offset + n_windows;
        segment_info.data{kF}.first_window_start_samp(kS) = seg_samples{kS}(1);
        total_window_offset = total_window_offset + n_windows;
    end
end

segment_info = explode_nested_tables(segment_info);

% Get top windows info for the requested channel
top_windows_all = inds_mfile.characteristic_sortinds;
top_windows_per_class = top_windows_all{day_ind}{chan_ind};
n_classes = length(top_windows_per_class);

mt_file_paths = cell(n_to_plot, n_classes);
lfp_inds = zeros(2, n_to_plot, n_classes);

figure;
tl = tiledlayout(n_classes, n_to_plot, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, sprintf('Example LFP per class (%s, %s)', rec_date, channel), 'Interpreter', 'none');

min_ymin = inf;
max_ymax = -inf;

% For each class:
for kK = 1:n_classes
    % Go through indices in order of "typicality" scores and find the first run within the
    % class that is long enough and includes the index without crossing a segment boundary
    top_windows = top_windows_per_class{kK};
    unexplored = true(length(top_windows), 1);
    start_wins = nan(n_to_plot, 1);
    n_found = 0;
    
    for kW = 1:length(top_windows)
        if ~unexplored(kW)
            continue;
        end
        
        this_win = top_windows(kW);
        kSeg = find(segment_info.last_window >= this_win, 1);
        
        if plot_windows == 1 % no need to use this complicated algorithm
            n_backward = 0;
            n_forward = 0;
        else
            seg_endpoint = segment_info.last_window(kSeg);
            forward_inds_to_look_for = [this_win + (1:ceil(plot_windows/2)-1), 0];
            n_forward = find(~ismember(forward_inds_to_look_for, top_windows), 1) - 1;
            n_forward = min(n_forward, seg_endpoint - this_win);

            seg_startpoint = segment_info.first_window(kSeg);
            backward_needed = plot_windows - n_forward - 1;
            backward_inds_to_look_for = [this_win - (1:backward_needed), 0];
            n_backward = find(~ismember(backward_inds_to_look_for, top_windows), 1) - 1;
            n_backward = min(n_backward, this_win - seg_startpoint);

            if n_backward < backward_needed && n_forward == length(forward_inds_to_look_for) - 1
                % haven't hit the edge yet, try pushing forward
                extra_forward_needed = plot_windows - (n_backward + n_forward + 1);
                extra_forward_inds = this_win + n_forward + (1:extra_forward_needed);
                if extra_forward_inds(end) <= seg_endpoint && all(ismember(extra_forward_inds, top_windows))
                    n_forward = n_forward + extra_forward_needed;
                end
            end
        end
            
        if n_backward + n_forward + 1 >= plot_windows
            n_found = n_found + 1;
            start_wins(n_found) = this_win - n_backward;
            if n_found == n_to_plot
                break;
            else
                % invalidate all windows that overlap, to avoid plotting the same data twice
                unexplored(top_windows >= this_win - n_backward - samps_in_window/2 & ...
                    top_windows <= this_win + n_forward + samps_in_window/2) = false;
            end
        else
            % invalidate just the current class run, to save time since it will never work
            unexplored(top_windows >= this_win - n_backward & top_windows <= this_win + n_forward) = false;
        end
    end
     
    if n_found < n_to_plot
        warning('Only %d of %d suitable LFP segments found for class %d', n_found, n_to_plot, kK);
    end
    
    % Now get info on where our sections of data are in terms of samples
    for kP = 1:n_to_plot
        nexttile;
        if kP > n_found
            continue;
        end
        
        start_win = start_wins(kP);
        kSeg = find(segment_info.last_window >= start_win, 1);
        kFile = segment_info.filenum(kSeg);
        mt_file_paths(kP, kK) = this_input_s.mt_res_in(kFile);
        window_offset = start_win - segment_info.first_window(kSeg);
        lfp_inds(1, kP, kK) = segment_info.first_window_start_samp(kSeg) + window_offset * step_samps;
        lfp_inds(2, kP, kK) = lfp_inds(1, kP, kK) - 1 + plot_samples;

        % Get LFP data corresponding to given time window and plot
        [~, rec_time] = fileparts(fileparts(mt_file_paths{kP, kK})); % sorry, hacky
        raw_path = fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s_%s.mat', rec_date, rec_time));
        raw_mfile = matfile(raw_path);
        time_inds = lfp_inds(1, kP, kK):lfp_inds(2, kP, kK);
        lfp_to_plot = organize_lfp(raw_mfile, chan_struct, [], [], [], time_inds);

        plot(1/Fs:1/Fs:plot_length, lfp_to_plot, 'k');
        this_ylims = get(gca, 'YLim');
        min_ymin = min(min_ymin, this_ylims(1));
        max_ymax = max(max_ymax, this_ylims(2));
        
        if kK == n_classes
            xlabel('Time (s)');
        end
        
        ylabel('\muV');
        
        if n_to_plot == 1
            title(sprintf('Class %d', kK));
        else
            title(sprintf('Class %d, sample %d', kK, kP));
        end
    end
end

for tilenum = 1:(n_classes * n_to_plot)
    nexttile(tilenum);
    ylim([min_ymin, max_ymax]);
end

end
