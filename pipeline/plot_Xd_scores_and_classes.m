function plot_Xd_scores_and_classes(nmf_mfile, X_chan_num)
% Given a NMF result file and a chosen channel whose space to focus on,
% plot the scores of each channel transformed to be comparable to that space,
% along with a class plot of the max scores after transformation.

if ischar(nmf_mfile)
    nmf_mfile = matfile(nmf_mfile);
end

% get relevant variables
nmf_U = nmf_mfile.nmf_U;
nmf_V = nmf_mfile.nmf_V;
kl_Xs = nmf_mfile.kl_Xs;
chan_names = nmf_mfile.chan_names;
n_chans = length(chan_names);
run_name = nmf_mfile.run_name;
time_axis = nmf_mfile.time_axis;

nmf_Vi = nmf_V{1};
nmf_Ui = nmf_U{1};
Xs_to_c = kl_Xs(:, X_chan_num);
ncomps = size(nmf_Vi{X_chan_num}, 2);

[~, Vs_i] = cellfun(@util.rescale_to_pmf, nmf_Ui, nmf_Vi, 'uni', false);
nmf_V_trans = cellfun(@(Vi, X) Vi * X, Vs_i, Xs_to_c, 'uni', false);

figure;
tl = tiledlayout(2, ceil(n_chans / 2), 'TileSpacing', 'compact');
title(tl, sprintf('Scores for %s transformed for min. KL divergence with %s', ...
    run_name, chan_names{X_chan_num}), 'Interpreter', 'none');

for kC = 1:length(chan_names)
    nexttile;
    sanePColor(time_axis, 1:ncomps, nmf_V_trans{kC}(:, 1:end-1).');
    xlabel('Time (s)');
    yticks(1:ncomps);
    ylabel('Component #');
    title(chan_names{kC}, 'Interpreter', 'none');
end

[~, maxclass_trans1] = cellfun(@(Vi) max(Vi(:, 1:end-1), [], 2), nmf_V_trans, 'uni', false);

figure;
class_plot(time_axis, chan_names, maxclass_trans1);
xlabel('Time (s)');
title('Max-component classes of transformed scores');

end
