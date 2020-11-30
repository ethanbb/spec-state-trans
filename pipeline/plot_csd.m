function plot_csd(recdate, probe_s, savedir)
% Plot and save CSD figure and data for selected date/probes.
% Inputs:
%   recdate: The date of interest, e.g. '2020-01-30'.
%            Uses all recordings with flashes on that date.
%   probe_s: Struct specifying which probes to use. Format:
%               Each key is a probe name (e.g. 'Probe1' if 'Probe1Name'
%                    is in the rec info)
%               Each value is the associated name/region, e.g. 'V1', 'M1'.
%   savedir: Where to save everything; defaults to fullfile(sr_dirs.results, recdate).
%
% Peaks of the resulting plots should be current *sinks*.

sr_dirs = prepSR;

% Get "snippit" files from this date
listing = dir(fullfile(sr_dirs.snippits, sprintf('snips_%s*', recdate)));
fnames = {listing.name};
rectimes = cellfun(@(fn) sscanf(fn, sprintf('snips_%s_%%[^.].mat', recdate)), ...
    fnames, 'uni', false);
files = fullfile({listing.folder}, fnames);
n_files = length(files);

probe_names = fieldnames(probe_s);
n_probes = length(probe_names);

snips = repmat({cell(n_files, 1)}, n_probes, 1);
probe_models = cell(n_probes, 1);

for kF = 1:n_files
    mfile = matfile(files{kF});
    if kF == 1
        Fs = mfile.finalSampR;
    else
        assert(mfile.finalSampR == Fs, 'Mismatched sample rates!');
    end
    
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
        
        snips{kP}{kF} = organize_lfp(mfile, struct(probe, 'all'));
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
            fhs(kP) = plot1csd(mean_snips{kP}{kF}, probe_s.(probe), probe_models{kP}, ...
                Fs, sprintf('%s_%s', recdate, rectimes{kF}));
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
    [fhs(kP), csd, time_axis, chan_axis, depth_axis] = ...
        plot1csd(mean_snips_combined{kP}, probe_s.(probe), probe_models{kP}, Fs, [recdate ' - all']);
    save(fullfile(savedir, ['csd_', probe_s.(probe), '.mat']), 'csd', 'time_axis', 'chan_axis', 'depth_axis');
end
savefig(fhs, fullfile(savedir, 'csd.fig'));
end

function [fh, gausCSD, time_axis, chan_axis, depth_axis] = ...
    plot1csd(chan_data, name, probe_model, Fs, title_note)

% get info about the probe
[num_chans, spacing] = util.get_probe_model_info(probe_model);
spacing = spacing / 1000; % to mm

% seconds of data per snippet before the stim
sec_pre = 1;

% setup for smoothed second spatial derivative
sigma = 2;
support = -ceil(3*sigma):ceil(3*sigma);

kernel = (-(sigma^2 - support.^2)/(sigma^5*sqrt(2*pi))) .* ...
    exp(-support.^2 / (2*sigma^2));

time = round(0.9*Fs:1.5*Fs); % 100 ms pre to 500 ms post
time_axis = time * 1000/Fs - 1000*sec_pre;

gausCSD = convn(chan_data, kernel.', 'valid');
chans_out = size(gausCSD, 1);

chans_removed = num_chans - chans_out;

fh = figure;
ax = axes;

chan_axis = chans_removed/2 + (1:chans_out);
depth_axis = chan_axis * spacing;

gausCSD = gausCSD(:, time);

sanePColor(time_axis, depth_axis, gausCSD);
axis tight;
shading interp;

ax.YDir = 'reverse';
ylim([depth_axis(1), depth_axis(end)]);
ylabel('Depth (mm)');

yyaxis right;
ax.YDir = 'reverse';
ylim([chan_axis(1), chan_axis(end)]);
ylabel('Channel # (on this probe)');

xlabel('Time since stimulus presentation (ms)');
titlestr = sprintf('CSD for %s', name);
if nargin > 3 && ~isempty(title_note)
    titlestr = [titlestr, ' (', title_note, ')'];
end
title(titlestr, 'Interpreter', 'none');

c = colorbar;
c.Label.String = 'CSD in mV/mm^2';
end
