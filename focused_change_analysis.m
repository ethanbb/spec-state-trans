%% get PCA-space change from file, or compute it

rng('shuffle');
rng_s = rng;

mfile = matfile('mt_res.mat');
regions = mfile.name;

if isprop(mfile, 'pca_change')
    smoothed_change = mfile.pca_change;
else
    smoothed_change = pca_change('mt_res.mat');
end

smoothed_change = num2cell(smoothed_change, 2);

%% find some points of high speed spectral change


lookdist = 3; % seconds back and forward to look for peak

types = {'real', 'sham'};

for kT = 1:2
    change_peaks.(types{kT}) = cell(1, 2);
end

bp = cell(1, 2);
prominence = cell(1, 2);
n_peaks = zeros(1, 2);

for kR = 1:2    
    change_bpeak = islocalmax(smoothed_change{kR}, ...
        'MinProminence', 0.9, ...
        'MinSeparation', lookdist * interp_factor);  % in seconds 

    change_peaks.real{kR} = find(change_bpeak);

    n_peaks(kR) = length(change_peaks.real{kR});
    
    change_peaks.sham{kR} = randperm(length(smoothed_change{kR}), n_peaks(kR));
end

%% do spike-triggered average of change in other region at these points

win_rel = -30*interp_factor:30*interp_factor; % 30 seconds before and after
win_rel_sec = win_rel / interp_factor;

blook = abs(win_rel_sec) <= lookdist;
win_rel_look_sec = win_rel_sec(blook);
    
for kT = 1:2
    type = types{kT};
    
    wins.(type) = cell(2, 2); % wins{i, j} = windows in region j based on peaks in i
    sta.(type) = cell(2, 2); % sta{i, j} = average of region j based on peaks in i
    peak_rel.(type) = nan(1, 2); % peak_rel(i) = average peak in other region aligned to i

    for kI = 1:2

        win_inds = change_peaks.(type){kI}(:) + win_rel; % creates npeaks x nwin matrix
        win_inds = max(win_inds, 1);
        win_inds = min(win_inds, length(smoothed_change{kI}));

        for kJ = 1:2
            wins.(type){kI, kJ} = smoothed_change{kJ}(win_inds);
            sta.(type){kI, kJ} = mean(wins.(type){kI, kJ});
        end

        % locate the average peak
%         kpeaks = find(islocalmax(sta.(type){kI, 3-kI}));
%         [~, kkpeak] = min(abs(win_rel(kpeaks))); % find peak closest to time 0
%         kpeak = kpeaks(kkpeak);
        % maybe just get the most prominent peak
        kpeak = find(islocalmax(sta.(type){kI, 3-kI}(blook), 'MaxNumExtrema', 1));
        if ~isempty(kpeak)            
            peak_rel.(type)(kI) = win_rel_look_sec(kpeak);
        end
    end
end

%% plot STA with peak

axs = gobjects(2, 2);
h_zero = gobjects(2, 2);
h_peak = gobjects(2, 2);
ylims = cell(1, 2);

sta_fig = figure;

for kR = 1:2
    
    axs(1, kR) = subplot(2, 2, kR);
    plot(win_rel_sec, sta.real{kR, kR}, 'k', 'LineWidth', 2);
    hold on;
    plot(win_rel_sec, sta.real{kR, 3-kR}, 'b', 'LineWidth', 2);
    ylims{kR} = get(gca, 'YLim');
    h_zero(1, kR) = plot([0, 0], ylims{kR}, 'k');
    h_peak(1, kR) = plot(peak_rel.real(kR) * [1, 1], ylims{kR}, 'b');
    axis tight;
    title(sprintf('%s spectral change aligned to %s spectral change peaks', ...
        regions{3-kR}, regions{kR}));
    
    axs(2, kR) = subplot(2, 2, kR + 2);
    plot(win_rel_sec, sta.sham{kR, kR}, 'k', 'LineWidth', 2);
    hold on;
    plot(win_rel_sec, sta.sham{kR, 3-kR}, 'b', 'LineWidth', 2);
    h_zero(2, kR) = plot([0, 0], ylims{kR}, 'k');
    h_peak(2, kR) = plot(peak_rel.sham(kR) * [1, 1], ylims{kR}, 'b');
    axis tight;
    title('Control - peak time selected at random');
end

%% bootstrap to get a distribution of peak lags

n_boot = 10000;
n_show = 20; % random ones to plot

hist_fig = figure;

for kT = 1:2
    type = types{kT};
    cis.(type) = cell(1, 2);
    
    for kR = 1:2
        k_show = randperm(n_boot, n_show);

        peak_rel_boot = nan(n_boot, 1);

        for kB = 1:n_boot
            samp = datasample(1:n_peaks(kR), n_peaks(kR));
            sta_boot = mean(wins.(type){kR, 3-kR}(samp, :));

            % locate average peak
%             kpeaks = find(islocalmax(sta_boot));
%             [~, kkpeak] = min(abs(win_rel(kpeaks))); % find peak closest to time 0
%             kpeak = kpeaks(kkpeak);
%             peak_rel_boot(kB) = win_rel_sec(kpeak);

            % record most prominent peak
            kpeak = find(islocalmax(sta_boot(blook), 'MaxNumExtrema', 1));
            if ~isempty(kpeak)
                peak_rel_boot(kB) = win_rel_look_sec(kpeak);
            end

            % plot a subset of bootstraps in dotted line
            if ismember(kB, k_show)
                plot(axs(kT, kR), win_rel_sec, sta_boot, 'b:');
            end
        end
        
        ylim(axs(kT, kR), ylims{kR});

        % plot confidence interval
        alpha = 0.95;
        ci_low = prctile(peak_rel_boot, 100*(1-alpha)/2);
        ci_high = prctile(peak_rel_boot, 100*(alpha + (1-alpha)/2));
        cis.(type){kR} = [ci_low, ci_high];

        h_ci = patch(axs(kT, kR), ...
            [ci_low * [1, 1], ci_high * [1, 1]], ...
            [ylims{kR}, ylims{kR}([2, 1])] + [-1, 1, 1, -1], ...
            'b', 'EdgeColor', 'b', 'LineStyle', '--', 'FaceAlpha', 0.1);

        if strcmp(type, 'real')
            legend(axs(kT, kR), [h_zero(kT, kR), h_peak(kT, kR), h_ci], ...
                sprintf('Peaks of %s change', regions{kR}), ...
                sprintf('Peak of mean aligned %s change', regions{3-kR}), ...
                sprintf('%.0f%% CI', alpha * 100));
        else
            legend(axs(kT, kR), [h_zero(kT, kR), h_peak(kT, kR), h_ci], ...
                sprintf('Random points in %s change', regions{kR}), ...
                sprintf('Peak of mean aligned %s change', regions{3-kR}), ...
                sprintf('%.0f%% CI', alpha * 100));
        end
        xlabel(axs(kT, kR), sprintf('Time relative to %s change peak (s)', regions{kR}));

        % plot histogram of bootstrapped peaks
        figure(hist_fig);
        subplot(2, 2, kR + 2*(kT-1));
        histogram(peak_rel_boot, 'Normalization', 'pdf');
        xlabel(sprintf('Time relative to %s peak (s)', regions{kR}));
        ylabel('Probability density');
        title(sprintf('Bootstrapped mean %s peak times (%s)', regions{3-kR}, type));
    end
end

%% Save

savefig([sta_fig, hist_fig], 'rel_change_times.fig');