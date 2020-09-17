function hf = plot_kldiv_mat(kl_div, chans, title_line2)
% Plot a grid of KL divergences between pairs of channels (as obtained by kl_divergence_analysis)
%
% Inputs:
%   kl_div: n_chans x n_chans matrix of NMF score KL divergences
%   chans: cell array of channel names
%   title_line2: If nonempty, print this below the general figure title.

mean_kl_div = mean(kl_div, 3);
hf = figure;
sanePColor(mean_kl_div);
set(gca, 'YDir', 'reverse');

title_line1 = 'NMF score KL divergence - $$\min_X \mathbf{E}[D_{KL}(V_j || V_iX)]$$';
if nargin < 3 || isempty(title_line2)
    title(title_line1, 'Interpreter', 'latex');
else
    title({title_line1, title_line2}, 'Interpreter', 'latex');
end

xticks(1:length(chans));
xticklabels(chans);
xtickangle(45);
yticks(1:length(chans));
yticklabels(chans);

set(gca, 'TickLabelInterpreter', 'none');

ylabel('Channel i');
xlabel('Channel j');

colorbar;

end
