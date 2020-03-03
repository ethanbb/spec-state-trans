% Using output from pca_vis, highlight points that are visually different
% (i.e. low power in low fequencies) on a 3d scatter plot.

prepSR;

date = '2020-02-06';
time = '16-01-00';
savedir = fullfile(results_dir, date, time);

mt_res = load(fullfile(savedir, 'mt_res.mat'));
data = mt_res.pxx_pca_1rec_vis{1};

% plot PCA data
hfig = figure;
hax = axes;
hax.Interactions = [zoomInteraction, rotateInteraction, rulerPanInteraction];

dscatter(data(1:3, :));
hold on;

% highlight manually identified points
hstarts = [1521, 2197, 3722, 4530, 5180, 6081, 7619];
hends = [1614, 2634, 3814, 4702, 5259, 6240, 7712];

Fs_mt = 10;
htimes = arrayfun(@(s, e) Fs_mt*s+1:Fs_mt*e, hstarts, hends, 'uni', false);
htimes = cell2mat(htimes);

scatter3(data(1, htimes), data(2, htimes), data(3, htimes), 20, 'r', 'o');

% save figure
savefig(hfig, fullfile(savedir, 'highlight_pca_pts.fig'), 'compact');