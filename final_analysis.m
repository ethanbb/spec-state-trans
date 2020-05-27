% Compile results from nnmf_clustering_all.m and look at conditional entropy and transition timing.

%% Set up variables

prepSR;

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
res_mfiles = cellfun(@(r) matfile(fullfile(results_dir, r, 'mt_res.mat')), recs, 'uni', false);

n_classes = 6;

%% Compile conditional entropy results

cond_ent_cell = cellfun(@(mf) mf.rank_cond_ent, res_mfiles, 'uni', false);
cond_ents = cat(3, cond_ent_cell{:});
med_cond_ent = median(cond_ents, 3);

chans = res_mfiles{1}.name;
n_chans = length(chans);

figure;
sanePColor(med_cond_ent);
set(gca, 'YDir', 'reverse');
xticks(1:n_chans);
xticklabels(chans);
xtickangle(45);
yticks(1:n_chans);
yticklabels(chans);
ylabel('R1');
xlabel('R2');
c = colorbar;
c.Label.String = 'Entropy (bits)';

title('Median conditional entropy of R2 given R1');

%% Repeat above for mutual information

mut_info_cell = cellfun(@(mf) mf.rank_mut_info, res_mfiles, 'uni', false);
mut_infos = cat(3, mut_info_cell{:});
med_mut_info = median(mut_infos, 3);

chans = res_mfiles{1}.name;
n_chans = length(chans);

figure;
sanePColor(med_mut_info);
set(gca, 'YDir', 'reverse');
xticks(1:n_chans);
xticklabels(chans);
xtickangle(45);
yticks(1:n_chans);
yticklabels(chans);
c = colorbar;
c.Label.String = 'Information (bits)';

title('Median mutual information of classes');


%% Collect mutual information entries of each type

types = {'Same channel', 'Near depth', 'Far depth', 'Cross-region'};
type_snames = cellfun(@matlab.lang.makeValidName, types, 'uni', false);

linear_inds = {};
linear_inds{1} = find(eye(8));
linear_inds{2} = linear_inds{1} + repmat([1; -1], 4, 1);
linear_inds{3} = sort([linear_inds{1}; linear_inds{2}]) + repmat([repmat(2, 4, 1); repmat(-2, 4, 1)], 2, 1);
linear_inds{4} = find((1:8)' > 4 & 1:8 <= 4 | (1:8)' <= 4 & 1:8 > 4);

minf_by_type = struct;
rec_vec = (0:n_recs-1) * 64;

for kT = 1:4
    minf_by_type.(type_snames{kT}) = squeeze(cellfun(@(mat) median(mat(linear_inds{kT})), num2cell(mut_infos, [1 2])));
    
%     type_inds = reshape(linear_inds{kT} + rec_vec, [], 1);
%     cent_by_type.(type_snames{kT}) = cond_ents(type_inds);
end

figure;
viol = violinplot(minf_by_type);
xticklabels(types);
ylabel('Mutual information (bits)');

colors = lines(5);
colors(3,:) = []; % skip 3rd since its yellow blends w/ parula colormap
for kT = 1:4
    viol(kT).ViolinColor = colors(kT, :);
end

hold on;
x_lims = get(gca, 'XLim');
plot(x_lims, log2(n_classes) * [1, 1], 'k--');
xlim(x_lims);
ylim([0, log2(n_classes)]);
curr_ticks = get(gca, 'YTick');
curr_labels = get(gca, 'YTickLabel');
yticks([curr_ticks(1:end-1), log2(n_classes)]);
yticklabels([curr_labels(1:end-1); {'max'}]);

title('Mutual information by pair type across recordings');

%% Do nonparametric stats (mutual info)

p1 = ranksum(minf_by_type.SameChannel, minf_by_type.NearDepth);
p2 = ranksum(minf_by_type.NearDepth, minf_by_type.FarDepth);
p3 = ranksum(minf_by_type.FarDepth, minf_by_type.Cross_region);

%% Conditional entropy vs. distance

elec_dist = 20; % in um

dists = zeros(numel(cond_ents) / 2 - n_recs * 8, 1);
cents = dists;
minfs = dists;

count = 0;
for kR = 1:n_recs
    res_opts = res_mfiles{kR}.options;
    
    for k1 = 1:8
        if k1 <= 4
            c2s = setdiff(1:4, k1);
        else
            c2s = setdiff(5:8, k1);
        end
        
        for k2 = c2s
            count = count + 1;
            dists(count) = elec_dist * abs(res_opts.chans(k1) - res_opts.chans(k2));
            cents(count) = cond_ents(k1, k2, kR);
            minfs(count) = mut_infos(k1, k2, kR);
        end
    end
end

assert(count == length(dists), 'Uh oh, counting error');

figure;
scatter(dists, cents, 'filled');
xlabel('Electrode distance (um)');
ylabel('Conditional entropy (bits)');

%% Unconditional entropy

entropies = zeros(n_recs, n_chans);
for kR = 1:n_recs
    classes = res_mfiles{kR}.rank_classes;
    for kC = 1:n_classes
        pC = sum(classes == kC) / size(classes, 1);
        new_entropies = entropies(kR, :) - pC .* log2(pC);
        
        % avoid nans
        bUpdate = pC > 0;        
        entropies(kR, bUpdate) = new_entropies(bUpdate);
    end
end

med_entropy = median(entropies);

figure;
bar(1:8, med_entropy);
xticklabels(chans);
xtickangle(45);
ylabel('Entropy (bits)');
title('Median entropy over recordings');

%% Transition timing comparison

% for each rec, make 3d cell of chan x from-class x to-class transition times.
all_trans = cell(n_recs, n_chans, n_classes, n_classes);
for kR = 1:n_recs
    times = res_mfiles{kR}.time_grid;
    full_classes = nan(length(times), n_chans);
    pxx1 = res_mfiles{kR}.pxx(1,1);
    b_nan = any(isnan(pxx1{1}));
    full_classes(~b_nan, :) = res_mfiles{kR}.rank_classes;
    
    for kC = 1:n_chans
        for kS = 2:size(full_classes, 1)
            k1 = full_classes(kS-1, kC);
            k2 = full_classes(kS, kC);
            
            if k1 ~= k2 && ~isnan(k1) && ~isnan(k2)
                all_trans{kR, kC, k1, k2}(end+1) = times(kS);
            end
        end
    end
end

%% Now use all_trans to quantify transition lags

max_dist = 60; % max dist b/w matching transitions in seconds

rec_matches = cell(n_classes, n_recs);
total_trans = 0;
total_matched = 0;

for kR = 1:n_recs
    chan_pair_matches = cell(n_classes, 16);
    for c1 = 1:4
        for c2 = 5:8
            
            chans = [c1, c2];
    
            trans = squeeze(all_trans(kR, :, :, :));
            n_classes = size(trans, 2);
            to_trans = cell(2, n_classes);
    
            % consolidate transition times *to* each class
            for kC = 1:2
                chan = chans(kC);
                for kK = 1:n_classes
                    to_trans{kC, kK} = sort(horzcat(trans{chan, setdiff(1:n_classes, kK), kK}));
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
            chan_pair_matches(:, (c1-1)*4 + (c2-4)) = k_matches;
        end
    end
    
    for kK = 1:n_classes
        rec_matches{kK, kR} = cell2mat(chan_pair_matches(kK,:));
    end
end
all_matches = horzcat(rec_matches{:});
match_diffs = -diff(all_matches);

figure;
t = tiledlayout(7, 1, 'TileSpacing', 'compact');
title(t, 'Lags of state transitions in MC from transitions to same state in V1');

h1 = nexttile([5, 1]);
histogram(match_diffs);
xlabel('Lag (s)');
ylabel('Count');
title(sprintf('%d/%d transitions matched (%.1f%%)', ...
    total_matched, total_trans, total_matched/total_trans*100));

h2 = nexttile([2, 1]);
boxplot(match_diffs, 'Orientation', 'horizontal', 'BoxStyle', 'filled', 'Colors', 'k');
box off;

linkaxes([h1, h2], 'x');

% Wilcoxon signed rank test:
%[p, h, stats] = signrank(match_diffs);

%% Look at differences across symmetric bins

figure;
lags = sort(match_diffs(match_diffs >= 0));
leads = sort(-match_diffs(match_diffs < 0));

h1 = plot(lags, 1:length(lags), 'b', 'LineWidth', 1.5);
hold on;
h2 = plot(leads, 1:length(leads), 'r', 'LineWidth', 1.5);

xlim([0, max(abs(match_diffs))]);
xlabel('Absolute timing difference (s)');
ylabel('Cumulative count');
legend([h1, h2], 'MC lags', 'MC leads', 'Location', 'northwest');

%% Bootstrap to get bounds on above figure

n_trans = length(match_diffs);
n_boot = 10000;
x_axis = 0:0.1:60;

lag_results = zeros(n_boot, length(x_axis));
lead_results = lag_results;

for kB = 1:n_boot
    sample = match_diffs(randi(n_trans, 1, n_trans));
    
    lags = sort(sample(sample >= 0));
    leads = sort(-sample(sample < 0));
    
    yvals_lag = 1:length(lags);
    % remove non-unique x points for interp1
    b_rem = lags == [lags(2:end), nan];
    lags = lags(~b_rem);
    yvals_lag = yvals_lag(~b_rem);
    
    yvals_lead = 1:length(leads);
    % remove non-unique x points for interp1
    b_rem = leads == [leads(2:end), nan];
    leads = leads(~b_rem);
    yvals_lead = yvals_lead(~b_rem);
    
    lag_results(kB, :) = interp1(lags, yvals_lag, x_axis);
    lead_results(kB, :) = interp1(leads, yvals_lead, x_axis);
end

hold on
h = area(x_axis, [min(lag_results)', (max(lag_results)-min(lag_results))']);
for kH = 1:2
    h(kH).EdgeColor = 'b';
    h(kH).LineStyle = '--';
    h(kH).FaceColor = 'b';
    h(kH).ShowBaseLine = false;
end
h(1).FaceAlpha = 0;
h(2).FaceAlpha = 0.3;

h = area(x_axis, [min(lead_results)', (max(lead_results)-min(lead_results))']);
for kH = 1:2
    h(kH).EdgeColor = 'r';
    h(kH).LineStyle = '--';
    h(kH).FaceColor = 'r';
    h(kH).ShowBaseLine = false;
end
h(1).FaceAlpha = 0;
h(2).FaceAlpha = 0.3;

% plot(x_axis, max(lag_results), 'b--');
% plot(x_axis, min(lag_results), 'b--');
% plot(x_axis, max(lead_results), 'r--');
% plot(x_axis, min(lead_results), 'r--');
legend([h1, h2], 'MC lags', 'MC leads', 'Location', 'northwest');
title('MC state transitions lagging or leading V1, by spread');
