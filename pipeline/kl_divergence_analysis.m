function kl_divergence_analysis(nmf_res_mfiles, channel_subsets)
% Compute directed transformations from NMF scores of each channel to each other channel within
% each day to minimize average KL divergence, then plot these KL divergences to look at
% near vs. distant vs. cross-region channel pairs.
%
% nmf_res_mfiles is a cell of writable MatFile objects,
% including null model data, on which to operate. 
% 
% channel_subsets, if not empty, should be a struct. Each value is a set of channels to operate on -
% either a vector of indices into nmf_V/nmf_U or a cell of channel names, which can include 
% wildcards (* and ?) to select multiple channels. (A channel may appear in
% more than one such set). The corresponding key identifies the set and will be prepended to the
% saved variables in each matfile. The default is to use all channels and the name 'kl'.

if exist('channel_subsets', 'var') && ~isempty(channel_subsets)
    assert(isstruct(channel_subsets) && isscalar(channel_subsets), ...
        'channel_subsets must be a scalar struct');
else
    channel_subsets = struct('kl', '*');
end

set_names = fieldnames(channel_subsets);
n_sets = length(set_names);

n_files = length(nmf_res_mfiles);
for kF = 1:n_files
    %% Get dataset info
    
    res_mfile = nmf_res_mfiles{kF};
    run_name = res_mfile.run_name;
    
    chans = res_mfile.chan_names;
        
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
    [~, Vs_i] = cellfun(@util.rescale_to_pmf, Us{1}, Vs_i, 'uni', false);
    [~, Vs_j, Vs_null] = cellfun(@util.rescale_to_pmf, Us{2}, Vs_j, Vs_null, 'uni', false);
    
    %% Loop over channel sets
    for kS = 1:n_sets
        %% setup loop
        
        set_name = set_names{kS};
        fprintf('Running set: %s\n', set_name);
        
        % Get actual channels depending on the format
        set_val = channel_subsets.(set_name);
        if ischar(set_val) || isstring(set_val)
            set_val = {set_val};
        end
        
        if isnumeric(set_val)
            set_chans = set_val;
        elseif iscell(set_val)
            % take channels based on names
            regexps = regexptranslate('wildcard', set_val);
            joined_regexp = strjoin(regexps, '|');
            
            set_chans = find(~cellfun('isempty', regexp(chans, joined_regexp)));
        else
            error('Unrecognized channel subset format');
        end
        
        n_chans = length(set_chans);
        
        Vs_i_set = Vs_i(set_chans);
        Vs_j_set = Vs_j(set_chans);
        Vs_null_set = Vs_null(set_chans);
        set_chan_names = chans(set_chans);        
        
        % element i,j = min_X(E[D_{KL}(Vj || Vi*X)]), where Vi is a random vector of NMF scores in channel i.
        kl_div = zeros(n_chans, n_chans, n_times);
        kl_div_null = zeros(n_chans, n_chans, n_times);
        Xs = cell(n_chans);
        Xs_null = cell(n_chans);
        
        % Keep KL divergence as vector over time for now
        kldiv_fn = @(Vi, Vj, X) sum(Vj .* log2((Vj+eps) ./ (Vi*X+eps)), 2);
        
        %% Find KL divergences        
        for iC = 1:n_chans
            Vi = Vs_i_set{iC};
            chan1_name = set_chan_names{iC};
            
            parfor jC = 1:n_chans                
                Vj = Vs_j_set{jC};
                
                % Status
                fprintf('\n---------- %s, %s -> %s, real ----------\n', run_name, chan1_name, set_chan_names{jC});
                
                Xs{iC, jC} = optimize_xentropy(Vi, Vj);
                
                % Save KL divergence
                kl_div(iC, jC, :) = kldiv_fn(Vi, Vj, Xs{iC, jC});
                
                % repeat for null model
                fprintf('\n---------- %s, %s -> %s, null ----------\n', run_name, chan1_name, set_chan_names{jC});
                Vj = Vs_null_set{jC};
                Xs_null{iC, jC} = optimize_xentropy(Vi, Vj);
                kl_div_null(iC, jC, :) = kldiv_fn(Vi, Vj, Xs_null{iC, jC});
            end
        end
        
        %% Save
        res_mfile.([set_name, '_Xs']) = Xs;
        res_mfile.([set_name, '_divs']) = kl_div;
        res_mfile.([set_name, '_Xs_null']) = Xs_null;
        res_mfile.([set_name, '_divs_null']) = kl_div_null;
        res_mfile.([set_name, '_chans']) = set_chan_names;
        
        %% Plot KL divergence
        mfile_dir = fileparts(res_mfile.Properties.Source);
        
        hf = plot_kldiv_mat(kl_div, set_chan_names, sprintf('%s - %s', run_name, set_name));
        savefig(hf, fullfile(mfile_dir, sprintf('score_DKL_%s_%s.fig', run_name, set_name)));
        
        hf2 = plot_kldiv_mat(kl_div_null, set_chan_names, sprintf('%s - %s (null model)', run_name, set_name));
        savefig(hf2, fullfile(mfile_dir, sprintf('score_DKL_%s_%s_null.fig', run_name, set_name)));
        
    end
end

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
