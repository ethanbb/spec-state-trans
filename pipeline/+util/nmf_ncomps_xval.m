function [n_comps, hfig] = nmf_ncomps_xval(data, min_comps, max_comps, pow_thresh)
% Do cross-validation to test how many NMF components to use.
% Inputs:
%   data:        the non-negative data matrix to decompose
%   min_comps:   smallest # of components to try
%   max_comps:   largest # of components to try
%   pow_thresh:  fraction of explained input power necessary to justify adding each component
% Outputs:
%   n_comps:     best # of components
%   hfig:        handle to figure showing cross-validation results

n_comp_options = min_comps:max_comps;
BCV_err = nnmf_k_ver(data, round(size(data) / 5), 5, min_comps, max_comps, [], 500000);

total_power = norm(data, 'fro');
explained = 1 - BCV_err / total_power;

% use all components that explain at least 1% of power
n_comps = find(diff(mean(explained)) < pow_thresh, 1);

if isempty(n_comps)
    n_comps = n_comp_options(end);
    warning('All components explain > 1% of power - consider trying more.');
end

if nargout > 1
    hfig = figure('Visible', 'off');
    scatter(reshape(repmat(n_comp_options, size(explained, 1), 1), [], 1), explained(:), 'k', 'filled');
    hold on;
    plot(n_comp_options, mean(explained), 'b', 'LineWidth', 1);
    xticks(n_comp_options);
    xlabel('Number of NMF components');
    ylabel('Fraction of explained power');
    title('NMF cross-validation');

    ylims = get(gca, 'YLim');
    h = plot([1 1] * (n_comps + 0.5), ylims, 'r--');
    ylim(ylims);
    legend(h, '1% power threshold', 'Location', 'southeast');
end

end
