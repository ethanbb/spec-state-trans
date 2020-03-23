% Compare extrema of spectral change velocity, consolidating data from all analyzed
% recordings (2020-01-30 - 2020-03-11).

recs = {
    '2020-01-30/16-03-00'
    '2020-01-31/12-52-00'
    '2020-01-31/15-26-00'
    '2020-02-06/13-47-00'
    '2020-02-06/16-01-00'
    '2020-03-05/12-50-00'
    '2020-03-05/14-50-00'
    '2020-03-06/12-55-00'
    '2020-03-06/14-56-00'
    '2020-03-10/12-57-00'
    '2020-03-10/14-19-00'
    '2020-03-11/12-31-00'
    '2020-03-11/14-32-00'
    };

n_recs = length(recs);

all_matched_peaks = cell(1, n_recs);
all_matched_troughs = cell(1, n_recs);

for kR = 1:n_recs
    res = matfile(fullfile(results_dir, recs{kR}, 'mt_res.mat'));
    all_matched_peaks{kR} = res.matched_peaks;
    all_matched_troughs{kR} = res.matched_troughs;
end

all_matched_peaks = cell2mat(all_matched_peaks);
all_matched_troughs = cell2mat(all_matched_troughs);

peak_diffs = -diff(all_matched_peaks);
trough_diffs = -diff(all_matched_troughs);

figure;
t = tiledlayout(7, 2, 'TileSpacing', 'compact');
title(t, 'Lags of spectral change extrema in MC from matched extrema in V1');

h1 = nexttile([5, 1]);
histogram(peak_diffs, 10, 'BinLimits', [-45, 45]);
xlabel('Lag (s)');
ylabel('Count');
title(sprintf('Change peaks (N = %d)', length(peak_diffs)));

h2 = nexttile([5, 1]);
histogram(trough_diffs, 10, 'BinLimits', [-45, 45]);
xlabel('Lag (s)');
ylabel('Count');
title(sprintf('Change troughs (N = %d)', length(trough_diffs)));

h3 = nexttile([2, 1]);
boxplot(peak_diffs, 'Orientation', 'horizontal', 'BoxStyle', 'filled', 'Colors', 'k');
box off;

h4 = nexttile([2, 1]);
boxplot(trough_diffs, 'Orientation', 'horizontal', 'BoxStyle', 'filled', 'Colors', 'k');
box off;

linkaxes([h1, h3], 'x');
linkaxes([h2, h4], 'x');

