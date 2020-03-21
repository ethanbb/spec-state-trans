function plot_csd(recdate)

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

for kF = 1:n_files
    fhs = gobjects(2, 1);
    
    for kR = 1:2
        gausCSD = convn(mean_snips{kF}(reg_inds{kR}, :), kernel.', 'valid');
        chans_out = size(gausCSD, 1);

        chans_removed = length(reg_inds{kR}) - chans_out;

        fhs(kR) = figure;

        chan_axis = (kR-1)*32 + chans_removed/2 + (1:chans_out);

        pcolor(time_axis, chan_axis, flipud(gausCSD(:, time)));
        shading interp;
        set(gca, 'YDir', 'reverse');
        axis tight;

        xlabel('Time since stimulus presentation (ms)');
        ylabel('Channel #');
        title(sprintf('CSD for %s (%s_%s)', regions{kR}, recdate, rectimes{kF}), ...
            'Interpreter', 'none');

        c = colorbar;
        c.Label.String = 'CSD in mV/mm^2';
    end
    
    savefig(fhs, fullfile(results_dir, recdate, sprintf('csd_%s.fig', rectimes{kF})));
end

end