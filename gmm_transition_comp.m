% Compare times of transitions to particular states between channels, consolidating data from all analyzed
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

max_dist = 60; % max dist b/w matching transitions in seconds

rec_matches = cell(1, n_recs);
total_trans = 0;
total_matched = 0;

for kR = 1:n_recs
    res = matfile(fullfile(results_dir, recs{kR}, 'mt_res.mat'));
    
    trans = res.gmm_transitions;
    n_classes = size(trans, 2);
    to_trans = cell(2, n_classes);
    
    % consolidate transition times *to* each class
    for kC = 1:2
       for kK = 1:n_classes
           to_trans{kC, kK} = sort(horzcat(trans{kC, setdiff(1:n_classes, kK), kK}));
       end
    end
    
    % now find matches
    k_matches = cell(1, n_classes);
    for kK = 1:n_classes
        nearest_trans = cell(2, 1);
        n_trans = cellfun('length', to_trans(:, kK));
        
        for kC = 1:2
            nearest_trans{kC} = zeros(1, n_trans(kC));
            for kT = 1:n_trans(kC)
                [dist, nearest] = min(abs(to_trans{kC, kK}(kT) - to_trans{3-kC, kK}));
                if ~isempty(nearest) && dist <= max_dist
                    nearest_trans{kC}(kT) = nearest;
                end
            end
        end
        
        % now that both classes have found their potential matches, loop again
        for kC = 1:2
            % set "nearest trans" entries for which the match is not reciprocal to 0
            % those that are already 0 will see the augmented 0 and thus will not match.
            augmented_matches = [0, nearest_trans{3-kC}];
            nearest_trans{kC}(augmented_matches(nearest_trans{kC} + 1) ~= 1:n_trans(kC)) = 0;
        end
        
        % add to totals
        n_matches = sum(nearest_trans{1} > 0);
        total_trans = total_trans + sum(n_trans);
        total_matched = total_matched + 2 * n_matches; % each match = 2 matched transitions
        
        k_matches{kK} = zeros(2, n_matches);
        for kC = 1:2
            k_matches{kK}(kC, :) = to_trans{kC, kK}(nearest_trans{3-kC}(nearest_trans{3-kC} > 0));
        end
    end
    rec_matches{kR} = cell2mat(k_matches);
end
all_matches = cell2mat(rec_matches);
match_diffs = -diff(all_matches);

figure;
t = tiledlayout(7, 2, 'TileSpacing', 'compact');
title(t, 'Lags of state transitions in MC from transitions to same state in V1');

h1 = nexttile([5, 1]);
histogram(match_diffs);
xlabel('Lag (s)');
ylabel('Count');
title(sprintf('%d/%d transitions matched (%.1f%%)', ...
    total_matched, total_trans, total_matched/total_trans*100));

h2 = nexttile([5, 1]);
qqplot(match_diffs);

h3 = nexttile([2, 1]);
boxplot(match_diffs, 'Orientation', 'horizontal', 'BoxStyle', 'filled', 'Colors', 'k');
box off;

linkaxes([h1, h3], 'x');

% statistics (reported in text box)
% 1-sample t-test:
[h, p, ci, stats] = ttest(match_diffs);
% power analysis:
mean_diff = mean(match_diffs);
std_diff = std(match_diffs);
nout = sampsizepwr('t', [0, std_diff], mean_diff, 0.9);
