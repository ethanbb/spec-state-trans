function kl_divergence_analysis(nmf_res_mfiles)
% Compute directed transformations from NMF scores of each channel to each other channel within
% each day to minimize average KL divergence, then plot these KL divergences to look at
% near vs. distant vs. cross-region channel pairs.
%
% nmf_res_mfiles is a cell of writable MatFile objects,
% including null model data, on which to operate. 

n_files = length(nmf_res_mfiles);
for kF = 1:n_files
    %% Get dataset info
    
    res_mfile = nmf_res_mfiles{kF};
    run_name = res_mfile.run_name;
    
    chans = res_mfile.chan_names;
    n_chans = numel(chans);
    
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
            fprintf('\n---------- %s, %s -> %s, real ----------\n', run_name, chans{iC}, chans{jC});
            
            Xs{iC, jC} = optimize_xentropy(Vi, Vj);
            
            % Save KL divergence
            kl_div(iC, jC, :) = kldiv_fn(Vi, Vj, Xs{iC, jC});
            
            % repeat for null model
            fprintf('\n---------- %s, %s -> %s, null ----------\n', run_name, chans{iC}, chans{jC});
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
    mfile_dir = fileparts(res_mfile.Properties.Source);
    
    hf = plot_kldiv_mat(kl_div, chans, run_name);
    savefig(hf, fullfile(mfile_dir, 'score_DKL.fig'));

    hf2 = plot_kldiv_mat(kl_div_null, chans, [run_name ' (null model)']);
    savefig(hf2, fullfile(mfile_dir, 'score_DKL_null.fig'));
end

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
