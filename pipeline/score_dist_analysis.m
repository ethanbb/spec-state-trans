function score_dist_analysis(nmf_res_mfiles, dist_type, channel_subsets)
% Compute directed transformations from NMF scores of each channel to each other channel within
% each day to minimize some type of distance between scores, then plot these distances to look at
% near vs. distant vs. cross-region channel pairs.
%
% nmf_res_mfiles is a cell of writable MatFile objects,
% including null model data, on which to operate. 
% 
% channel_subsets, if not empty, should be a struct. Each value is a set of channels to operate on -
% either a vector of indices into nmf_V/nmf_U or a cell of channel names, which can include 
% wildcards (* and ?) to select multiple channels. (A channel may appear in
% more than one such set). The corresponding key identifies the set and will be prepended to the
% saved variables in each matfile, along with the dist_type.
% The default is to use all channels and the name 'all'.

arguments
    nmf_res_mfiles (:,1) cell
    dist_type (1,:) char {mustBeMember(dist_type, {'kl_div', 'L2_dist', 'recon_err'})} = 'kl_div'
    channel_subsets (1,1) struct = struct('all', '*')
end

% change things that depend on the distance type
recon_err_fn = @(Vi, Vj, Uj, X) norm(Vj*Uj' - Vi*X*Uj', 'fro');

switch dist_type
    case 'kl_div'
        make_prob_dist = true;
        optim_fn = @(Vi, Vj, Uj, X) -sum(sum(Vj .* log2(Vi*X+eps)));
        result_fn = @(Vi, Vj, Uj, X) mean(sum(Vj .* log2((Vj+eps) ./ (Vi*X+eps)), 2));
        
    case 'L2_dist'
        make_prob_dist = false;
        optim_fn = @(Vi, Vj, Uj, X) norm(Vj - Vi*X, 'fro');
        result_fn = optim_fn;
        
    case 'recon_err' % minimizing reconstruction error directly - empirically seems quite close to euclidean
        make_prob_dist = false;
        optim_fn = recon_err_fn;
        result_fn = recon_err_fn;
        
    otherwise
        error('Unrecognized distance type');
end

set_names = fieldnames(channel_subsets);
n_sets = length(set_names);

n_files = length(nmf_res_mfiles);
for kF = 1:n_files
    %% Get dataset info
    
    res_mfile = nmf_res_mfiles{kF};
    run_name = res_mfile.run_name;
    
    chans = res_mfile.chan_names;
        
    Vs = res_mfile.nmf_V;
    Us = res_mfile.nmf_U;
    
    % Split into 2 repetitions - one used as Q/i and one as P/j in KL divergence calculation.
    % This allows a fair lower bound when comparing each channel with itself.
    Us_i = Us{1};
    Us_j = Us{2};
    Vs_i = Vs{1};
    Vs_j = Vs{2};
    Vs_null = res_mfile.null_V;
    Vs_null = Vs_null{1}; % Only use on one side for now

    % rescale all to be probability distributions
    if make_prob_dist
        [~, Vs_i, Vs_null] = cellfun(@util.rescale_to_pmf, Us_i, Vs_i, Vs_null, 'uni', false);
        [Us_j, Vs_j] = cellfun(@util.rescale_to_pmf, Us_j, Vs_j, 'uni', false);
    end
    
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

        Us_j_set = Us_j(set_chans);
        Vs_i_set = Vs_i(set_chans);
        Vs_j_set = Vs_j(set_chans);
        Vs_null_set = Vs_null(set_chans);
        set_chan_names = chans(set_chans);      
        
        % element i,j = min_X(E[D_{KL}(Vj || Vi*X)]), where Vi is a random vector of NMF scores in channel i.
        dist = zeros(n_chans, n_chans);
        dist_null = zeros(n_chans, n_chans);
        recon_err = zeros(n_chans, n_chans);
        recon_err_null = zeros(n_chans, n_chans);
        Xs = cell(n_chans);
        Xs_null = cell(n_chans);
        
        %% Find KL divergences        
        for iC = 1:n_chans
            Vi_real = Vs_i_set{iC};
            Vi_null = Vs_null_set{iC};
            chan1_name = set_chan_names{iC};
            
            parfor jC = 1:n_chans
                Uj = Us_j_set{jC};
                Vj = Vs_j_set{jC};
                
                % Status
                fprintf('\n---------- %s, %s -> %s, real ----------\n', run_name, chan1_name, set_chan_names{jC});
                
                Xs{iC, jC} = optimize_dist_fn(Vi_real, Vj, Uj, optim_fn, make_prob_dist);
                
                % Save KL divergence and reconstruction error
                dist(iC, jC) = result_fn(Vi_real, Vj, Uj, Xs{iC, jC}); %#ok<PFBNS>
                recon_err(iC, jC) = recon_err_fn(Vi_real, Vj, Uj, Xs{iC, jC});
                
                % repeat for null model
                fprintf('\n---------- %s, %s -> %s, null ----------\n', run_name, chan1_name, set_chan_names{jC});
                Xs_null{iC, jC} = optimize_dist_fn(Vi_null, Vj, Uj, optim_fn, make_prob_dist);
                dist_null(iC, jC) = result_fn(Vi_null, Vj, Uj, Xs_null{iC, jC});
                recon_err_null(iC, jC) = recon_err_fn(Vi_null, Vj, Uj, Xs{iC, jC});
            end
        end
        
        %% Save
        var_prefix = [set_name, dist_type];
        res_mfile.([var_prefix, '_Xs']) = Xs;
        res_mfile.([var_prefix, '_Xs_null']) = Xs_null;
        res_mfile.([var_prefix, 's']) = dist;
        res_mfile.([var_prefix, 's_null']) = dist_null;
        res_mfile.([var_prefix, '_recon_err']) = recon_err;
        res_mfile.([var_prefix, '_recon_err_null']) = recon_err_null;
        
        res_mfile.([set_name, '_chans']) = set_chan_names;
        
        %% Plot KL divergence
        mfile_dir = fileparts(res_mfile.Properties.Source);
        
        hf = plot_dist_mat(dist, set_chan_names, ...
            sprintf('%s - %s', run_name, set_name), dist_type);
        savefig(hf, fullfile(mfile_dir, ...
            sprintf('score_%s_%s_%s.fig', dist_type, run_name, set_name)));
        
        hf2 = plot_dist_mat(dist_null, set_chan_names, ...
            sprintf('%s - %s (null model)', run_name, set_name), dist_type);
        savefig(hf2, fullfile(mfile_dir, ...
            sprintf('score_%s_%s_%s_null.fig', dist_type, run_name, set_name)));
        
        if ~strcmp(dist_type, 'recon_err')
            recon_type = ['recon_from_', dist_type];
            
            % Also plot reconstruction error
            hf3 = plot_dist_mat(recon_err, set_chan_names, ...
                sprintf('%s - %s', run_name, set_name), recon_type);
            savefig(hf3, fullfile(mfile_dir, ...
                sprintf('score_%s_%s_%s.fig', recon_type, run_name, set_name)));
            
            hf4 = plot_dist_mat(recon_err_null, set_chan_names, ...
                sprintf('%s - %s (null model)', run_name, set_name), recon_type);
            savefig(hf4, fullfile(mfile_dir, ...
                sprintf('score_%s_%s_%s_null.fig', recon_type, run_name, set_name)));
        end
        
    end
end

end

function X = optimize_dist_fn(Vi, Vj, Uj, optim_fn, use_prob_dist)
% Get the matrix X that minimizes cross-entropy of Vi*X relative to Vj.

ki = size(Vi, 2);
kj = size(Vj, 2);
kmin = min(ki, kj);

% fmincon problem parameters
fn = @(X) optim_fn(Vi, Vj, Uj, X);
x0 = ones(ki, kj) / kj;
x0(1:kmin, 1:kmin) = 0.8*eye(kmin)*kmin/kj + 0.2*ones(kmin)/kj; % a bit arbitrary
if use_prob_dist
    Aeq = repmat(eye(ki), 1, kj); % restricts each row to sum to 1
    beq = ones(ki, 1);
else
    Aeq = [];
    beq = [];
end
lb = zeros(size(x0));
ub = ones(size(x0));
optimopts = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e5);

X = fmincon(fn, x0, [], [], Aeq, beq, lb, ub, [], optimopts);

end
