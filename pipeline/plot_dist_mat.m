function hf = plot_dist_mat(dist_mats, chan_names, title_line2, dist_type, plot_type, chan_names_with_regions)
% Plot a grid of score distances between pairs of channels (as obtained by score_dist_analysis)
%
% Inputs:
%   kl_div: n_chans x n_chans matrix of NMF score KL divergences
%   chan_names: cell array of channel names
%   title_line2: If nonempty, print this below the general figure title.
%   dist_type: what kind of distance this is - see switch statement below for options
%   plot_type (optional): one of:
%       'full' (default): plot the full matrix
%       'lower': lower triangle, including diagonal
%       'lower_nodiag': lower triangle, not including diagonal
%       'upper', 'upper_nodiag'
%   chan_names_with_regions (optional): if provided, used to determine what regions channels are in
%                                       and make text labels accordingly. Regions should be the 
%                                       text of each channel name up to the first underscore.

n_chans = length(chan_names);
assert(all(size(dist_mats, 1:2) == n_chans), 'Mismatch in number of channels/names');

mean_dist_mat = mean(dist_mats, 3, 'omitnan');

% nan-out entries depending on plot type
if ~exist('plot_type', 'var') || isempty(plot_type)
    plot_type = 'full';
end

cols = 1:n_chans;
rows = cols';
switch plot_type
    case 'full' % do nothing
    case 'lower'
        mean_dist_mat(cols > rows) = nan;
    case 'lower_nodiag'
        mean_dist_mat(cols >= rows) = nan;
        rows = rows(2:end);
        cols = cols(1:end-1);
    case 'upper'
        mean_dist_mat(cols < rows) = nan;
    case 'upper_nodiag'
        mean_dist_mat(cols <= rows) = nan;
        rows = rows(1:end-1);
        cols = cols(2:end);
    otherwise
        error('Unrecognized plot type');
end

hf = figure('Position', [300, 300, 600, 500]);
sanePColor(mean_dist_mat);
hold on;
xlim([cols(1)-0.5, cols(end)+0.5]);
ylim([rows(1)-0.5, rows(end)+0.5]);
set(gca, 'YDir', 'reverse');

if ~exist('dist_type', 'var') || isempty(dist_type)
    dist_type = 'kl_div';
end

interpreter = 'latex';
chan1_label = 'Channel 1 (depth in um)';
chan2_label = 'Channel 2 (depth in um)';

switch dist_type
    case 'kl_div'
        title_line1 = 'NMF score KL divergence - $$\min_X \mathbf{E}[D_{KL}(V_j || V_iX)]$$';
        chan1_label = 'Channel i (depth in $$\mu$$m)';
        chan2_label = 'Channel j (depth in $$\mu$$m)';
    case 'L2_dist'
        title_line1 = 'Euclidean score distance - $$\min_X \|V_j - V_iX\|_F$$';
        chan1_label = 'Channel i (depth in $$\mu$$m)';
        chan2_label = 'Channel j (depth in $$\mu$$m)';
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
    case 'norm_mutual_info'
        title_line1 = 'Norm. mut. info. of classes between channels';
        interpreter = 'none';
    case 'norm_mutual_info_z'
        title_line1= 'NMI of classes between channels - z-score against null model';
        interpreter = 'none';
    case 'trans_sync_scores'
        title_line1 = 'Discrete class transition mean synchrony';
        interpreter = 'none';
    case 'trans_sync_scores_z'
        title_line1 = 'Discrete class transition synchrony - z-score against null model';
        interpreter = 'none';
    case 'cca'
        title_line1 = 'Mean canonical correlation of NMF scores';
        interpreter = 'none';
    case 'cca_z'
        title_line1 = 'Canonical correlation of NMF scores - z-score against null model';
        interpreter = 'none';
    otherwise
        error('Unrecognized matrix type');
end

if strcmp(interpreter, 'latex')
    title_line2 = strrep(title_line2, '_', '\_');
end

if ~exist('title_line2', 'var') || isempty(title_line2)
    title(title_line1, 'Interpreter', interpreter);
else
    title({title_line1, title_line2}, 'Interpreter', interpreter);
end

chan_names = strrep(chan_names, '_', '\_');

xticks(cols);
xticklabels(chan_names(cols));
xtickangle(45);
yticks(rows);
yticklabels(chan_names(rows));
box off;

yl = ylabel(chan1_label);
xl = xlabel(chan2_label);

colorbar;

if exist('chan_names_with_regions', 'var') && ~isempty(chan_names_with_regions)
    % add region labels  
    regions = string(strtok(chan_names_with_regions, '_'));
    
    ylabel(sprintf('%s\n\n', chan1_label));
    xlabel(sprintf('\n%s', chan2_label));
    
    ax_pos = get(gca, 'Position');
    
    ytext_pos = yl.Extent(1) + yl.Extent(3);
    yl.Units = 'normalized';
    ytext_pos_norm = ax_pos(1) + (yl.Extent(1)+yl.Extent(3)) * ax_pos(3);
    make_ytext = @(loc, t) text(ytext_pos, loc, t, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none', ...
        'HorizontalAlignment', 'right');
    
    xtext_pos = xl.Extent(2) - xl.Extent(4) + 0.2;
    xl.Units = 'normalized';
    xtext_pos_norm = ax_pos(2) + (xl.Extent(2)+xl.Extent(4)) * ax_pos(4);
    make_xtext = @(loc, t) text(loc, xtext_pos, t, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');
    
    base_pos_x = cols(1) - 0.5;
    base_pos_y = rows(1) - 0.5;
    curr_reg = regions(1);
    
    data_to_y_norm = @(d) ax_pos(2) + (length(rows)-(d-rows(1)+0.5))*ax_pos(4)/length(rows);
    data_to_x_norm = @(d) ax_pos(1) + (d-cols(1)+0.5)*ax_pos(3)/length(cols);

    for kC = 2:n_chans + 1
        % do y and x axes separately, in case the diagonal isn't included
        if kC > rows(1) && kC <= rows(end)+1 && (kC == rows(end)+1 || regions(kC) ~= curr_reg)
            end_pos = kC - 0.5;
            mid_pos = (base_pos_y + end_pos) / 2;
            make_ytext(mid_pos, curr_reg);
            
            if kC <= rows(end)
                % dotted line between tick labels
                ypos_norm = data_to_y_norm(end_pos);
                annotation('line', ytext_pos_norm - [0.05, 0], [ypos_norm, ypos_norm], ...
                    'LineStyle', ':', 'LineWidth', 1);
                
                % heavy line separating regions in plot
                switch plot_type
                    case 'full'
                        sep_x = [cols(1)-0.5, cols(end)+0.5];
                    case {'upper', 'upper_nodiag'}
                        sep_x = [end_pos, cols(end)+0.5];
                    case {'lower', 'lower_nodiag'}
                        sep_x = [cols(1)-0.5, end_pos];
                end
                plot(sep_x, [end_pos, end_pos], 'k', 'LineWidth', 1.5);
                
                base_pos_y = end_pos;
            end
        end
        
        if kC > cols(1) && kC <= cols(end)+1 && (kC == cols(end)+1 || regions(kC) ~= curr_reg)
            end_pos = kC - 0.5;
            mid_pos = (base_pos_x + end_pos) / 2;
            make_xtext(mid_pos, curr_reg);
            
            if kC <= cols(end)
                % dotted line between tick labels
                xpos_norm = data_to_x_norm(end_pos);
                annotation('line', [xpos_norm, xpos_norm], xtext_pos_norm + [0.01, -0.04], ...
                    'LineStyle', ':', 'LineWidth', 1);
                
                % heavy line separating regions in plot
                switch plot_type
                    case 'full'
                        sep_y = [rows(1)-0.5, rows(end)+0.5];
                    case {'upper', 'upper_nodiag'}
                        sep_y = [rows(1)-0.5, end_pos];
                    case {'lower', 'lower_nodiag'}
                        sep_y = [end_pos, rows(end)+0.5];
                end
                plot([end_pos, end_pos], sep_y, 'k', 'LineWidth', 1.5);
                
                base_pos_x = end_pos;
            end
        end
        
        if kC <= n_chans && regions(kC) ~= curr_reg
            curr_reg = regions(kC);        
        end
    end
end

end
