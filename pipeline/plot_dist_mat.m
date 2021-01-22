function hf = plot_dist_mat(dist_mats, chans, title_line2, dist_type)
% Plot a grid of KL divergences between pairs of channels (as obtained by kl_divergence_analysis)
%
% Inputs:
%   kl_div: n_chans x n_chans matrix of NMF score KL divergences
%   chans: cell array of channel names
%   title_line2: If nonempty, print this below the general figure title.

mean_dist_mat = mean(dist_mats, 3);
hf = figure;
sanePColor(mean_dist_mat);
set(gca, 'YDir', 'reverse');

if nargin < 4 || isempty(dist_type)
    dist_type = 'kl_div';
end

switch dist_type
    case 'kl_div'
        title_line1 = 'NMF score KL divergence - $$\min_X \mathbf{E}[D_{KL}(V_j || V_iX)]$$';
    case 'L2_dist'
        title_line1 = 'Norm of score difference - $$\min_X \|V_j - V_iX\|_F$$';
    case 'recon_err'
        title_line1 = 'Spectrogram reconstruction error of channel j from channel i';
    case 'recon_from_kl_div'
        title_line1 = 'Reconstruction error from minimizing score KL divergence';
    case 'recon_from_L2_dist'
        title_line1 = 'Reconstruction error from minimizing norm of score difference';
end

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
