function hf = plot_dist_mat(dist_mats, chans, title_line2, dist_type)
% Plot a grid of score distances between pairs of channels (as obtained by score_dist_analysis)
%
% Inputs:
%   kl_div: n_chans x n_chans matrix of NMF score KL divergences
%   chans: cell array of channel names
%   title_line2: If nonempty, print this below the general figure title.
%   dist_type: what kind of distance this is - see switch statement below for options

mean_dist_mat = mean(dist_mats, 3);
hf = figure;
sanePColor(mean_dist_mat);
set(gca, 'YDir', 'reverse');

if nargin < 4 || isempty(dist_type)
    dist_type = 'kl_div';
end

interpreter = 'latex';
cmap = 'parula';
switch dist_type
    case 'kl_div'
        title_line1 = 'NMF score KL divergence - $$\min_X \mathbf{E}[D_{KL}(V_j || V_iX)]$$';
    case 'L2_dist'
        title_line1 = 'Euclidean score distance - $$\min_X \|V_j - V_iX\|_F$$';
    case 'L2_dist_from_null'
        title_line1 = 'Euclidean score proximity compared to null model';
    case 'recon_err'
        title_line1 = 'Spectrogram reconstruction error of channel j from channel i';
        interpreter = 'none';
    case 'recon_from_kl_div'
        title_line1 = 'Reconstruction error from minimizing score KL divergence';
        interpreter = 'none';
    case 'recon_from_L2_dist'
        title_line1 = 'Reconstruction error from minimizing norm of score difference';
        interpreter = 'none';
    case 'mutual_info'
        title_line1 = 'Mutual information of classes between channels';
        interpreter = 'none';
        cmap = 'jet';
    case 'norm_mutual_info'
        title_line1 = 'Norm. mut. info. of classes between channels';
        interpreter = 'none';
        cmap = 'jet';
    otherwise
        error('Unrecognized matrix type');
end

if strcmp(interpreter, 'latex')
    title_line2 = strrep(title_line2, '_', '\_');
end

if nargin < 3 || isempty(title_line2)
    title(title_line1, 'Interpreter', interpreter);
else
    title({title_line1, title_line2}, 'Interpreter', interpreter);
end

chans = strrep(chans, '_', '\_');

xticks(1:length(chans));
xticklabels(chans);
xtickangle(45);
yticks(1:length(chans));
yticklabels(chans);

ylabel('Channel i');
xlabel('Channel j');

colormap(cmap);
colorbar;

end
