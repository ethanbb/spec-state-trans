function typicality_sortinds = sort_class_data_by_typicality(nmf_res_path_or_mfile, interactive)
% Given NMF results for a single rat (day), produces an array of one cell per channel;
% these cells each contain an array of one cell per class;
% the inner cells contain a permutation of the indices into the concatenated data that belong to
% that class, ordered from closest to the class centroid to farthest away, by Mahalanobis distance.
% If 'interactive' is true, pause to show plot of clusters for each channel.

if ~exist('interactive', 'var') || isempty(interactive)
    interactive = false;
end

if ischar(nmf_res_path_or_mfile) || isStringScalar(nmf_res_path_or_mfile)
    nmf_mfile = matfile(nmf_res_path_or_mfile);
else
    nmf_mfile = nmf_res_path_or_mfile;
    nmf_mfile.Properties.Writable = false;
end

% Get classes and (normalized, concatenated) data
n_classes = nmf_mfile.nmf_comps;
n_chans = length(n_classes);
classes = nmf_mfile.filtered_classes;
classes = classes{1};  % n_chans x 1 cell
data = nmf_mfile.pxx_cat;
chan_names = nmf_mfile.chan_names;

typicality_sortinds = cell(n_chans, 1);
for kC = 1:n_chans
    this_n_classes = n_classes(kC);
    this_data_centered = data{kC} - mean(data{kC}, 2);
    typicality_sortinds{kC} = cell(this_n_classes, 1);
    
    if interactive
        fig = figure;
        legend('Location', 'northeastoutside');
        hold on;
        [pca_loadings, pca_scores] = pca(this_data_centered', 'NumComponents', 2);
        pc_xbounds = [min(pca_scores(:, 1)), max(pca_scores(:, 1))];
        pc_ybounds = [min(pca_scores(:, 2)), max(pca_scores(:, 2))];
    end
    
    for kK = 1:this_n_classes
        this_class_inds = find(classes{kC} == kK);
        this_class_data = this_data_centered(:, this_class_inds)';
        mahal_dists = mahal(this_class_data, this_class_data);
        [~, mahal_sortinds] = sort(mahal_dists);
        typicality_sortinds{kC}{kK} = this_class_inds(mahal_sortinds);
        
        if interactive
            this_class_proj = this_class_data * pca_loadings;
            plot(this_class_proj(:, 1), this_class_proj(:, 2), '.', 'DisplayName', sprintf('Class %d', kK));
            mu = mean(this_class_data);
            sigma = cov(this_class_data);
            gausspdf = @(x, y) reshape(mvnpdf([x(:), y(:)] * pca_loadings', mu, sigma), size(x));
            fcontour(gausspdf, [pc_xbounds, pc_ybounds], 'DisplayName', sprintf('Class %d Gaussian', kK));
        end
    end
    
    if interactive
        xlabel('PC 1');
        ylabel('PC 2');
        title(sprintf('%s classes (close to continue)', chan_names{kC}), 'Interpreter', 'none');
        waitfor(fig);
    end
end

end
