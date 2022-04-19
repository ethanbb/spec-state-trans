function [typicality_sortinds, n_most_typical_spaced] = sort_class_data_by_typicality(...
    nmf_res_path_or_mfile, raw_path_or_mfile, use_autocorr, interactive)
% Given NMF results for a single rat (day), produces an array of one cell per channel;
% these cells each contain an array of one cell per class;
% the inner cells contain a permutation of the indices into the concatenated data that belong to
% that class, ordered from most to least "distinctive" of that class (determined by ratio of the
% class's NMF score to sum of scores of all classes). 
% Second output has a similar structure, but gives only up to 100 indices for each class,
% and ensures that they are spaced far enough apart in time (within each class) to avoid
% autocorrelation (based on the raw data).
% If 'interactive' is true, pause to show plot of clusters for each channel.

if ~exist('use_autocorr', 'var') || isempty(use_autocorr)
    use_autocorr = false;
end

if ~exist('interactive', 'var') || isempty(interactive)
    interactive = true;
end

if ischar(nmf_res_path_or_mfile) || isStringScalar(nmf_res_path_or_mfile)
    nmf_mfile = matfile(nmf_res_path_or_mfile);
else
    nmf_mfile = nmf_res_path_or_mfile;
    nmf_mfile.Properties.Writable = false;
end

if use_autocorr
    if ischar(raw_path_or_mfile) || isStringScalar(raw_path_or_mfile)
        raw_mfile = matfile(raw_path_or_mfile);
    else
        raw_mfile = raw_path_or_mfile;
        raw_mfile.Properties.Writable = false;
    end

    % get autocorr distance max over a few channels
    chans_to_use = 8:16:120;
    raw_data = raw_mfile.meanSubFullTrace(chans_to_use, :);
    autocorr_lengths = zeros(length(chans_to_use), 1);
    start_nlags = 2000;
    for kC = 1:length(chans_to_use)
        nlags = start_nlags;
        while true
            [acf, lags, bounds] = autocorr(raw_data(kC, :), 'NumLags', nlags);
            first_null = find(acf < bounds(1) & acf > bounds(2), 1);
            if ~isempty(first_null)
                autocorr_lengths(kC) = lags(first_null) / raw_mfile.finalSampR;
                break;
            else
                nlags = nlags * 2;
            end
        end
    end
    autocorr_length = max(autocorr_lengths);
    autocorr_mt_samps = ceil(autocorr_length * 10); % step size = 0.1
else
    autocorr_mt_samps = 0;
end

% Get classes and (normalized, concatenated) data
n_classes = nmf_mfile.nmf_comps;
n_chans = length(n_classes);
classes = nmf_mfile.filtered_classes;
classes = classes{1};  % n_chans x 1 cell
nmf_scores = nmf_mfile.nmf_V;
nmf_scores = nmf_scores{1}; % n_chans x 1 cell
data = nmf_mfile.pxx_cat;
chan_names = nmf_mfile.chan_names;

typicality_sortinds = cell(n_chans, 1);
n_most_typical_spaced = cell(n_chans, 1);
for kC = 1:n_chans   
    this_n_classes = n_classes(kC);
    this_data_centered = data{kC} - mean(data{kC}, 2);
    typicality_sortinds{kC} = cell(this_n_classes, 1);
    n_most_typical_spaced{kC} = cell(this_n_classes, 1);
    
    % calculate typicality measure based on NMF scores
    this_scores = nmf_scores{kC}; % time x modes
    scores_rel = this_scores ./ sum(this_scores, 2);
    
    if interactive
        fig = figure;
        nax = 2;
        axs = gobjects(nax, 1);
        for kax = 1:nax
            axs(kax) = subplot(nax, 1, kax);
            legend('Location', 'northeastoutside');
            xlabel('PC 1');
            ylabel('PC 2');
            zlabel('PC 3');
            box on;
            hold on;
        end

        ncomp = 3;
        [~, pca_scores] = pca(this_data_centered', 'NumComponents', ncomp);
%         [~, pca_scores] = pca(scores_rel, 'NumComponents', ncomp);  % for 3D visualization of relative scores
%         pc_bounds = [min(pca_scores(:, 1:ncomp))', max(pca_scores(:, 1:ncomp))'];
    end
    
    for kK = 1:this_n_classes
        this_class_inds = find(classes{kC} == kK);
%         this_class_data = this_data_centered(:, this_class_inds)';
%         mahal_dists = mahal(this_class_data, this_class_data);
%         [~, mahal_sortinds] = sort(mahal_dists);
%         typicality_sortinds{kC}{kK} = this_class_inds(mahal_sortinds);

        this_class_scores_rel = scores_rel(this_class_inds, kK);
        [~, score_rel_sortinds] = sort(this_class_scores_rel, 'descend');
        abs_sortinds = this_class_inds(score_rel_sortinds);
        
        % Get the top 100, with minimum spacing based on autocorrelation
        n_to_take = 100;
        
        if use_autocorr
            valid_to_take = true(length(this_class_inds), 1);
            most_typical_inds = zeros(n_to_take, 1);
            n_taken = 0;
            for k = 1:length(abs_sortinds)
                if ~valid_to_take(k)
                    continue;
                end

                n_taken = n_taken + 1;
                this_ind = abs_sortinds(k);
                most_typical_inds(n_taken) = this_ind;
                valid_to_take(abs_sortinds >= this_ind - autocorr_mt_samps & abs_sortinds <= this_ind + autocorr_mt_samps) = false;
                if n_taken == n_to_take
                    break;
                end
            end

            if n_taken < n_to_take
                warning('%s, class %d: only found %d of %d spaced-apart points', chan_names{kC}, kK, n_taken, n_to_take);
                most_typical_inds = most_typical_inds(1:n_taken);
            end
        else  % just take the first 100 points
            if length(abs_sortinds) < n_to_take
                warning('%s, class %d: only found %d of %d points', chan_names{kC}, kK, length(abs_sortinds), n_to_take);
                n_to_take = length(abs_sortinds);
            end
            
            most_typical_inds = abs_sortinds(1:n_to_take);
            k = n_to_take;
        end
        
        typicality_sortinds{kC}{kK} = abs_sortinds;
        n_most_typical_spaced{kC}{kK} = most_typical_inds;        
        
        if interactive
            this_class_pca_scores = pca_scores(this_class_inds, :);
            % plot(this_class_scores(:, 1), this_class_scores(:, 2), '.', 'DisplayName', sprintf('Class %d', kK));
            plot3(axs(1), this_class_pca_scores(:, 1), this_class_pca_scores(:, 2), this_class_pca_scores(:, 3), '.', 'DisplayName', sprintf('Class %d', kK));
%             mu = mean(this_class_data);
%             sigma = cov(this_class_data);
%             gausspdf = @(x, y) reshape(mvnpdf([x(:), y(:)] * pca_loadings', mu, sigma), size(x));
%             fcontour(gausspdf, [pc_xbounds, pc_ybounds], 'DisplayName', sprintf('Class %d Gaussian', kK));
            
            plot3(axs(2), pca_scores(most_typical_inds, 1), ...
                pca_scores(most_typical_inds, 2), ...
                pca_scores(most_typical_inds, 3), '.', 'DisplayName', sprintf('Class %d', kK));

            figure;
            histogram(this_class_scores_rel);
            hold on;
            thresh = this_class_scores_rel(score_rel_sortinds(k));
            plot([thresh, thresh], get(gca, 'YLim'), 'r--');
            title(sprintf('Class %d relative score distribution', kK));
        end
    end
    
    if interactive
        title(axs(1), sprintf('PCA of %s class spectra (close to continue)', chan_names{kC}), 'Interpreter', 'none');
        title(axs(2), sprintf('PCA of %s class spectra (%d most typical per class)', chan_names{kC}, n_to_take), 'Interpreter', 'none');
        waitfor(fig);
    end
end
end
