% Compute directed transformations from NMF scores of each channel to each other channel within
% each day to minimize average KL divergence, then plot these KL divergences to look at
% near vs. distant vs. cross-region channel pairs.

%% Set up variables

prepSR;

days = {
    '2020-01-30'
    '2020-01-31'
    '2020-02-06'
    '2020-03-05'
    '2020-03-06'
    '2020-03-10'
    '2020-03-11'
    };

n_days = length(days);

%% Loop through recordings
for kD = 1:n_days
    %% Get dataset info
    day = days{kD};

    res_path = fullfile(results_dir, day, 'nmf_res.mat');
    res_mfile = matfile(res_path, 'Writable', true);
    
    chans = res_mfile.chan_names;
    chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false);
    n_chans = numel(chans);
    
    freq_axis = res_mfile.freq_axis;
    time_axis = res_mfile.time_axis;
    n_times = length(time_axis);

    Vs = res_mfile.nmf_V;
    Us = res_mfile.nmf_U;

    % Split into 2 repetitions - one used as Q/i and one as P/j in KL divergence calculation.
    % This allows a fair lower bound when comparing each channel with itself.
    Vs_i = Vs{1};
    Vs_j = Vs{2};
    Vs_null = res_mfile.null_V;
    Vs_null = Vs_null{2}; % Only use on one side for now

    % rescale all to be probability distributions
    [~, Vs_i] = cellfun(@rescale_to_pmf, Us{1}, Vs_i, 'uni', false);
    [~, Vs_j, Vs_null] = cellfun(@rescale_to_pmf, Us{2}, Vs_j, Vs_null, 'uni', false);
    
    %% setup loop
    
    % element i,j = min_X(E[D_{KL}(Vj || Vi*X)]), where Vi is a random vector of NMF scores in channel i.
    kl_div = zeros(n_chans, n_chans, n_times);
    kl_div_null = zeros(n_chans, n_chans, n_times);
    Xs = cell(n_chans);
    Xs_null = cell(n_chans);
    
    % Keep KL divergence as vector over time for now
    kldiv_fn = @(Vi, Vj, X) sum(Vj .* log2((Vj+eps) ./ (Vi*X+eps)), 2);
        
    %% Find KL divergences
    for iC = 1:n_chans
        for jC = 1:n_chans
            %%
            Vi = Vs_i{iC};
            Vj = Vs_j{jC};
            
            % Status
            fprintf('\n---------- %s, %s -> %s, real ----------\n', day, chans{iC}, chans{jC});
            
            Xs{iC, jC} = optimize_xentropy(Vi, Vj);
            
            % Save KL divergence
            kl_div(iC, jC, :) = kldiv_fn(Vi, Vj, Xs{iC, jC});
            
            % repeat for null model
            fprintf('\n---------- %s, %s -> %s, null ----------\n', day, chans{iC}, chans{jC});
            Vj = Vs_null{jC};
            Xs_null{iC, jC} = optimize_xentropy(Vi, Vj);
            kl_div_null(iC, jC, :) = kldiv_fn(Vi, Vj, Xs_null{iC, jC});
        end
    end
    
    %% Save
    res_mfile.kl_Xs = Xs;
    res_mfile.kl_divs = kl_div;
    res_mfile.kl_Xs_null = Xs_null;
    res_mfile.kl_divs_null = kl_div_null;


    %% Plot KL divergence
    hf = plot_kldiv_mat(kl_div, chans, day);
    savefig(hf, fullfile(results_dir, day, 'score_DKL.fig'));

    hf2 = plot_kldiv_mat(kl_div_null, chans, [day ' (null model)']);
    savefig(hf2, fullfile(results_dir, day, 'score_DKL_null.fig'));
end

%% Get and plot median KL divergence for each channel pair (here still assuming V1 super...etc. labels are comparable)
% For each day, use only chan names specified in the csd results file (since I went back and
% filtered out channels that don't actually correspond to the layer they're supposed to be)

layers = {'L2/3', 'L4', 'L5'};
n_chans = length(layers) * 2;

all_mean_kl_divs = nan(n_chans, n_chans, n_days);
for kD = 1:n_days
    res_kld = load(fullfile(results_dir, days{kD}, 'nmf_res.mat'), 'kl_divs', 'chan_names');
    csd_chans_V1_s = load(fullfile(results_dir, days{kD}, 'csd_V1.mat'), 'chan_names');
    csd_chans_MC_s = load(fullfile(results_dir, days{kD}, 'csd_MC.mat'), 'chan_names');

    b_chan = [ismember(layers, csd_chans_V1_s.chan_names), ...
              ismember(layers, csd_chans_MC_s.chan_names)];

    all_mean_kl_divs(b_chan, b_chan, kD) = mean(res_kld.kl_divs(b_chan, b_chan, :), 3);
end

% little hacky, grab chan names again here in case I didn't run the loop above
chans = res_kld.chan_names;

med_kl_divs = nanmedian(all_mean_kl_divs, 3);
hf = plot_kldiv_mat(med_kl_divs, chans, sprintf('Median over %d days', n_days));
savefig(hf, fullfile(results_dir, 'res_figs', 'med_kl_div_days.fig'));

%% Try same thing but with difference from null model

all_mean_kl_divs_null = zeros(n_chans, n_chans, n_days);
for kD = 1:n_days
    res_mfile = matfile(fullfile(results_dir, days{kD}, 'nmf_res.mat'));
    all_mean_kl_divs_null(:, :, kD) = mean(res_mfile.kl_divs_null, 3);
end

med_kl_div_from_null = nanmedian(all_mean_kl_divs_null - all_mean_kl_divs, 3);
hf = plot_kldiv_mat(med_kl_div_from_null, chans, 'Amount below divergence from null model (median)');
savefig(hf, fullfile(results_dir, 'res_figs', 'med_kl_div_fromnull_days.fig'));

%% Make violin plot

types = {'Same channel', 'Same region', 'Cross-region'};
type_snames = cellfun(@matlab.lang.makeValidName, types, 'uni', false);
n_types = length(types);

linear_inds = cell(n_types, 1);
linear_inds{1} = find(eye(n_chans));
linear_inds{2} = find(blkdiag(ones(n_chans/2), ones(n_chans/2)) - eye(n_chans));
linear_inds{end} = setdiff(1:n_chans^2, vertcat(linear_inds{1:end-1}));

kl_divs_bytype = struct;
kl_divs_null_bytype = struct;

for kT = 1:n_types
    %mean_kldiv_oftype = @(mat) mean(mat(linear_inds{kT}));

    kl_divs_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_mean_kl_divs, [1, 2]), 'uni', false);
    kl_divs_bytype.(type_snames{kT}) = vertcat(kl_divs_perday{:});

    kl_divs_null_perday = cellfun(@(mat) mat(linear_inds{kT}), num2cell(all_mean_kl_divs_null, [1, 2]), 'uni', false);
    kl_divs_null_bytype.(type_snames{kT}) = vertcat(kl_divs_null_perday{:});
end

hf = figure; hold on;
viol1 = violinplot(kl_divs_bytype, [], 'ViolinColor', [1, 0, 0]);
viol2 = violinplot(kl_divs_null_bytype, [], 'ViolinColor', [0.5, 0.5, 0.5]);
legend([viol1(1).ViolinPlot, viol2(1).ViolinPlot], 'Real data', 'Markov-chain generated', 'Location', 'southeast');
xticklabels(types);
ylabel('KL divergence (bits)');

title({'Mean KL divergence of aligned NMF scores between pairs of channels,', ...
    'with second channel either real data or resampled based on discrete Markov chain'});

savefig(hf, fullfile(results_dir, 'res_figs', 'kl_div_withnull_violin.fig'));

%%
function hf = plot_kldiv_mat(kl_div, chans, title_line2)

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

function [new_U, varargout] = rescale_to_pmf(U, varargin)
% Makes V or cell of Vs (from NMF) into a series of probability mass functions by the following
% steps:
%  * Divide each V by the maximum sum over components over all Vs (such that they each sum to <= 1
%    at each point), and multiply U by the same factor.
%  * Add a column to U consisting of all zeros
%  * Add a column to each V consisting of 1 - sum(V, 2) (such that all rows sum to 1).

scale_factor = max(cellfun(@(v) max(sum(v, 2)), varargin));
scale_factor = scale_factor + eps(scale_factor); % make sure there's no discrepancy when dividing each element
Vs = cellfun(@(v) v / scale_factor, varargin, 'uni', false);
U = U * scale_factor;

new_U = [U, zeros(size(U, 1), 1)];
varargout = cellfun(@(v) [v, 1 - sum(v, 2)], Vs, 'uni', false);

end

function X = optimize_xentropy(Vi, Vj)
% Get the matrix X that minimizes cross-entropy of Vi*X relative to Vj.

ki = size(Vi, 2);
kj = size(Vj, 2);
kmin = min(ki, kj);

% fmincon problem parameters
fn = @(X) -sum(sum(Vj .* log2(Vi*X+eps)));
x0 = ones(ki, kj) / kj;
x0(1:kmin, 1:kmin) = 0.8*eye(kmin)*kmin/kj + 0.2*ones(kmin)/kj; % a bit arbitrary
Aeq = repmat(eye(ki), 1, kj); % restricts each row to sum to 1
beq = ones(ki, 1);
lb = zeros(size(x0));
ub = ones(size(x0));
optimopts = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e5);

X = fmincon(fn, x0, [], [], Aeq, beq, lb, ub, [], optimopts);

end