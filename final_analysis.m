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

%% Collect conditional entropy entries of each type

types = {'Same channel', 'Near depth', 'Far depth', 'Cross-region'};
type_snames = cellfun(@matlab.lang.makeValidName, types, 'uni', false);

linear_inds = {};
linear_inds{1} = find(eye(8));
linear_inds{2} = linear_inds{1} + repmat([1; -1], 4, 1);
linear_inds{3} = sort([linear_inds{1}; linear_inds{2}]) + repmat([repmat(2, 4, 1); repmat(-2, 4, 1)], 2, 1);
linear_inds{4} = find((1:8)' > 4 & 1:8 <= 4 | (1:8)' <= 4 & 1:8 > 4);

cent_by_type = struct;
rec_vec = (0:n_recs-1) * 64;

for kT = 1:4
    cent_by_type.(type_snames{kT}) = squeeze(cellfun(@(mat) median(mat(linear_inds{kT})), num2cell(cond_ents, [1 2])));
    
%     type_inds = reshape(linear_inds{kT} + rec_vec, [], 1);
%     cent_by_type.(type_snames{kT}) = cond_ents(type_inds);
end

figure;
viol = violinplot(cent_by_type);
xticklabels(types);
ylabel('Conditional entropy (bits)');

colors = lines(4);
for kT = 1:4
    viol(kT).ViolinColor = colors(kT, :);
end

hold on;
x_lims = get(gca, 'XLim');
plot(x_lims, log2(6) * [1, 1], 'k--');
xlim(x_lims);
ylim([0, log2(6)]);
curr_ticks = get(gca, 'YTick');
curr_labels = get(gca, 'YTickLabel');
yticks([curr_ticks(1:end-1), log2(6)]);
yticklabels([curr_labels(1:end-1); {'max'}]);

title('Conditional entropy by pair type across recordings');

%% Do nonparametric stats

p1 = ranksum(cent_by_type.SameChannel, cent_by_type.NearDepth);
p2 = ranksum(cent_by_type.NearDepth, cent_by_type.FarDepth);
p3 = ranksum(cent_by_type.FarDepth, cent_by_type.Cross_region);