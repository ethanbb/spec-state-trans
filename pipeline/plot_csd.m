function plot_csd(recdate, probe_s, bad_chan_s, savedir, artifact_table, stim_event_chan_s)
% Plot and save CSD figure and data for selected date/probes.
% Inputs:
%   recdate: The date of interest, e.g. '2020-01-30'.
%            Uses all recordings with flashes on that date.
%   probe_s: Struct specifying which probes to use. Format:
%               Each key is a probe name (e.g. 'Probe1' if 'Probe1Name'
%                    is in the rec info)
%               Each value is the associated name/region, e.g. 'V1', 'M1'.
%   bad_chans: A struct with the same keys as probe_s, with each value a vector of
%              channel numbers to exclude from the CSD (defaulting to none).
%   savedir: Where to save everything; defaults to fullfile(sr_dirs.results, recdate).
%   artifact_table: If not empty, table containing n_artifacts x 2 matrix of artifact times to
%                   exclude (in s) in variable 'artifacts', with recording names as the row names.
%   stim_event_chan_s: Event channel index for stimulation for each probe (used to determine
%                      snippit times for artifact removal) (defaults to channel 1)
%
% Peaks of the resulting plots should be current *sinks*.

sr_dirs = prepSR;

% Get "snippit" files from this date
listing = dir(fullfile(sr_dirs.snippits, sprintf('snips_%s*', recdate)));
fnames = {listing.name};
recnames = cellfun(@(fn) sscanf(fn, 'snips_%[^.].mat'), fnames, 'uni', false);
recname_parts = split(recnames(:), '_');
rectimes = recname_parts(:, 2);
files = fullfile({listing.folder}, fnames);
n_files = length(files);
assert(n_files > 0, ['No snippits found on ', recdate]);

probe_names = fieldnames(probe_s);
n_probes = length(probe_names);

if ~exist('bad_chan_s', 'var') || isempty(bad_chan_s)
    bad_chan_s = struct;
    for kP = 1:n_probes
        bad_chan_s.(probe_names{kP}) = [];
    end
end

if ~exist('stim_event_chan_s', 'var') || isempty(stim_event_chan_s)
    stim_event_chan_s = struct;
    for kP = 1:n_probes
        stim_event_chan_s.(probe_names{kP}) = 1;
    end
end

snips = repmat({cell(n_files, 1)}, n_probes, 1);
probe_models = cell(n_probes, 1);

for kF = 1:n_files
    mfile = matfile(files{kF});
    if kF == 1
        Fs = mfile.finalSampR;
    else
        assert(mfile.finalSampR == Fs, 'Mismatched sample rates!');
    end
    
    all_start_times = mfile.allStartTimes;
    
    for kP = 1:n_probes
        probe = probe_names{kP};
        
        % get model info
        ifo = mfile.info;
        model = ifo.([probe, 'Name']);
        if kF == 1
            probe_models{kP} = model;
        else
            assert(strcmp(model, probe_models{kP}), 'Mismatched probe models!');
        end
        
        snips{kP}{kF} = organize_lfp(mfile, struct(probe, 'all'), [], [], true);
        
        if exist('artifact_table', 'var')
            % exclude artifacts
            event_chan = stim_event_chan_s.(probe);
            start_times = all_start_times{event_chan};
            n_snippits = size(snips{kP}{kF}, 2);
            artifacts = artifact_table.artifacts{recnames{kF}};
            b_include_snippits = util.get_snippits_to_include(start_times, n_snippits, 1, 3, artifacts);
            snips{kP}{kF}(:, ~b_include_snippits, :) = [];
        end
    end
end

% take average for each recording
mean_snips = cellfun(@(probe_snips) ...
    cellfun(@(s) squeeze(mean(s, 2)), probe_snips, 'uni', false), ...
    snips, 'uni', false); % chans x samples

% make sure there's somewhere to save them
if ~exist('savedir', 'var') || isempty(savedir)
    savedir = fullfile(sr_dirs.results, recdate);
end

if ~exist(savedir, 'dir')
    if ~mkdir(savedir)
        error('Could not create folder %s to save results', savedir);
    end
end

% Do CSD for each file individually
if n_files > 1
    for kF = 1:n_files
        fhs = gobjects(n_probes, 1);

        for kP = 1:n_probes
            probe = probe_names{kP};
            snips_masked = mean_snips{kP}{kF};
            snips_masked(bad_chan_s.(probe), :) = nan;
            fhs(kP) = plot1csd(snips_masked, probe_s.(probe), probe_models{kP}, Fs, recnames{kF});
        end

        savefig(fhs, fullfile(savedir, sprintf('csd_%s.fig', rectimes{kF})));
    end
end

% Also do a CSD for all data combined
snips_combined = cellfun(@(probe_snips) cat(2, probe_snips{:}), snips, 'uni', false);
mean_snips_combined = cellfun(@(probe_snips) squeeze(mean(probe_snips, 2)), ...
    snips_combined, 'uni', false);

fhs = gobjects(n_probes, 1);
for kP = 1:n_probes
    probe = probe_names{kP};
    snips_masked = mean_snips_combined{kP};
    snips_masked(bad_chan_s.(probe), :) = nan;
    [fhs(kP), csd, time_axis, chan_axis, depth_axis] = ...
        plot1csd(snips_masked, probe_s.(probe), probe_models{kP}, Fs, [recdate ' - all']);
    save(fullfile(savedir, ['csd_', probe_s.(probe), '.mat']), 'csd', 'time_axis', 'chan_axis', 'depth_axis');
end
savefig(fhs, fullfile(savedir, 'csd.fig'));
end

function [fh, gausCSD, time_axis, chan_axis, depth_axis] = ...
    plot1csd(chan_data, name, probe_model, Fs, title_note)

% get info about the probe
[num_chans, spacing] = util.get_probe_model_info(probe_model);

% seconds of data per snippet before the stim
sec_pre = 1;

% setup for smoothed second spatial derivative
sigma = 2;
support = -ceil(3*sigma):ceil(3*sigma);
n_edge = (length(support)-1) / 2;
chan_axis = 1+n_edge:num_chans-n_edge;
depth_axis = chan_axis * spacing / 1000; % (in mm)

kernel = (-(sigma^2 - support.^2)/(sigma^5*sqrt(2*pi))) .* ...
    exp(-support.^2 / (2*sigma^2));

time = round(0.9*Fs:1.5*Fs); % 100 ms pre to 500 ms post
time_axis = time * 1000/Fs - 1000*sec_pre;
chan_data = chan_data(:, time);

% do the 2nd spatial derivative
gausCSD = nanconv(chan_data, kernel.', 'nanout');
gausCSD = gausCSD(chan_axis, :);  % "valid" channels

% do a version for display that interpolates over nans
bad_chans = find(any(isnan(chan_data), 2));
good_chans = setdiff(1:num_chans, bad_chans);

if ~isempty(good_chans)
    [X, Y] = meshgrid(time, good_chans);
    [Xq, Yq] = meshgrid(time, bad_chans);
    chan_data_interp = chan_data;
    chan_data_interp(bad_chans, :) = interp2(X, Y, chan_data(good_chans, :), Xq, Yq);

    gausCSD_interp = convn(chan_data_interp, kernel.', 'valid');
else
    % no good channels! (e.g. all artifacts)
    gausCSD_interp = gausCSD;
end

fh = figure;
ax = axes;

sanePColor(time_axis, chan_axis, gausCSD_interp);
axis tight;
shading interp;

ax.YDir = 'reverse';
ylim([chan_axis(1), chan_axis(end)]);
ylabel(sprintf('Channel # (%g {\\mu}m per channel)', spacing));

xlabel('Time since stimulus presentation (ms)');
titlestr = sprintf('CSD for %s', name);
if nargin > 3 && ~isempty(title_note)
    titlestr = [titlestr, ' (', title_note, ')'];
end
title(titlestr, 'Interpreter', 'none');

c = colorbar;
c.Label.String = 'CSD in mV/mm^2';

end
