function plot_csd(recdate)
% Peaks of the resulting plots should be current *sinks*.

prepSR;

% Get "snippit" files from this date
listing = dir(fullfile(snippits_dir, sprintf('snips_%s*', recdate)));
fnames = {listing.name};
rectimes = cellfun(@(fn) sscanf(fn, sprintf('snips_%s_%%[^.].mat', recdate)), ...
    fnames, 'uni', false);
files = fullfile({listing.folder}, fnames);
n_files = length(files);

snips = cell(n_files, 1);

for kF = 1:n_files
    temp = load(files{kF}, 'meanSubData', 'finalSampR', 'info');
    snips{kF} = organize_lfp(temp); % now 1:32 are visual, 33:64 are motor
    
    if kF == 1
        Fs = temp.finalSampR;
    else
        assert(temp.finalSampR == Fs, 'Mismatched sample rates!');
    end
end

% take average for each recording
mean_snips = cellfun(@(s) squeeze(mean(s, 2)), snips, 'uni', false); % chans x samples

% make sure there's somewhere to save them
savedir = fullfile(results_dir, recdate);
if ~exist(savedir, 'dir')
    % try to create it

    if ~mkdir(savedir)
        error('Could not create folder %s to save results', savedir);
    end
end

% Do CSD for each file individually
if n_files > 1
    for kF = 1:n_files
        fhs = gobjects(2, 1);

        for kR = 1:2
            fhs(kR) = plot1csd(mean_snips{kF}, kR, Fs, sprintf('%s_%s', recdate, rectimes{kF}));
        end

        savefig(fhs, fullfile(savedir, sprintf('csd_%s.fig', rectimes{kF})));
    end
end

% Also do a CSD for all data combined
snips_combined = cat(2, snips{:});
mean_snips_combined = squeeze(mean(snips_combined, 2));

fhs = gobjects(2, 1);
regions = {'V1', 'MC'};
for kR = 1:2
    [fhs(kR), csd, time_axis, chan_axis, depth_axis] = plot1csd(mean_snips_combined, kR, Fs, [recdate ' - all']);
    save(fullfile(savedir, ['csd_', regions{kR}, '.mat']), 'csd', 'time_axis', 'chan_axis', 'depth_axis');
end
savefig(fhs, fullfile(savedir, 'csd.fig'));
end

function [fh, gausCSD, time_axis, chan_axis, depth_axis] = plot1csd(chan_data, kR, Fs, title_note)
% Do CSD and plot
regions = {'V1', 'MC'};
reg_inds = {1:32, 33:64};
sec_pre = 1;

% setup for smoothed second spatial derivative
sigma = 2;
support = -ceil(3*sigma):ceil(3*sigma);

kernel = (-(sigma^2 - support.^2)/(sigma^5*sqrt(2*pi))) .* ...
    exp(-support.^2 / (2*sigma^2));

time = round(0.9*Fs:1.5*Fs); % 100 ms pre to 500 ms post
time_axis = time * 1000/Fs - 1000*sec_pre;

spacing = 25/1000; % in mm

gausCSD = convn(chan_data(reg_inds{kR}, :), kernel.', 'valid');
chans_out = size(gausCSD, 1);

chans_removed = length(reg_inds{kR}) - chans_out;

fh = figure;
ax = axes;

chan_axis = (kR-1)*32 + chans_removed/2 + (1:chans_out);
depth_axis = (chan_axis - (kR-1)*32) * spacing;

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
ylabel('Channel #');

xlabel('Time since stimulus presentation (ms)');
ylabel('Channel #');
titlestr = sprintf('CSD for %s', regions{kR});
if nargin > 3 && ~isempty(title_note)
    titlestr = [titlestr, ' (', title_note, ')'];
end
title(titlestr, 'Interpreter', 'none');

c = colorbar;
c.Label.String = 'CSD in mV/mm^2';
end