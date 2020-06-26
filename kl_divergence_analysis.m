% Compute directed transformations from NMF scores of each channel to each other channel within
% a recording to minimize average KL divergence, then plot these KL divergences to look at
% near vs. distant vs. cross-region channel pairs.

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

%% Loop through recordings
for kR = 1:n_recs
    %% Get dataset info
    rec = recs{kR};

    res_path = fullfile(results_dir, rec, 'mt_res.mat');
    res_mfile = matfile(res_path, 'Writable', true);
    
    mt_opts = res_mfile.options; % necessary due to MatFile quirk
    Fw = 1 / mt_opts.winstep; % "window rate"
    chan_inds = mt_opts.chans;
    
    chans = res_mfile.name;
    chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false);
    n_chans = numel(chans);
    
    freq_axis = res_mfile.freq_grid;
    time_axis = res_mfile.time_grid;
    
    %% setup loop
    Vs = res_mfile.pmf_V;
    
    % element i,j = min_X(E[D_{KL}(Vj || Vi*X)]), where Vi is a random vector of NMF scores in channel i.
    kl_div = zeros(n_chans);
    Xs = cell(n_chans);
    
    cross_ent_fn = @(Vi, Vj, X) -sum(sum(Vj .* log2(Vi*X+eps))); % function to minimize
    kl_div_fn = @(Vi, Vj, X) mean(sum(Vj .* log2((Vj+eps) ./ (Vi*X+eps)), 2)); % value to save
    
    % set diagonal Xs to identity
    for kC = 1:n_chans
        Xs{kC, kC} = eye(size(Vs{kC}, 2));
    end
        
    %% Find KL divergences
    for iC = 1:n_chans
        for jC = setdiff(1:n_chans, iC)
            %%
            Vi = Vs{iC};
            Vj = Vs{jC};
            
            ki = size(Vi, 2);
            kj = size(Vj, 2);
            kmin = min(ki, kj);
            
            % fmincon problem parameters
            fn = @(X) cross_ent_fn(Vi, Vj, X);
            x0 = ones(ki, kj) / kj;
            x0(1:kmin, 1:kmin) = 0.8*eye(kmin)*kmin/kj + 0.2*ones(kmin)/kj; % a bit arbitrary
            Aeq = repmat(eye(ki), 1, kj); % restricts each row to sum to 1
            beq = ones(ki, 1);
            lb = zeros(size(x0));
            ub = ones(size(x0));
            optimopts = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e5);
            
            Xs{iC, jC} = fmincon(fn, x0, [], [], Aeq, beq, lb, ub, [], optimopts);
            
            kl_div(iC, jC) = kl_div_fn(Vi, Vj, Xs{iC, jC});
        end
    end
    
    %% Save
    res_mfile.kl_Xs = Xs;
    res_mfile.kl_divs = kl_div;
    
    %% Plot KL divergence
    hf = plot_kldiv_mat(kl_div, chans, rec);    
    savefig(hf, fullfile(results_dir, rec, 'score_DKL.fig'));
end

%% Get and plot median KL divergence for each channel pair (here still assuming V1 super...etc. labels are comparable)

all_kl_divs = zeros(n_chans, n_chans, n_recs);
for kR = 1:n_recs
    res_mfile = matfile(fullfile(results_dir, recs{kR}, 'mt_res.mat'));
    all_kl_divs(:, :, kR) = res_mfile.kl_divs;
end

med_kl_divs = median(all_kl_divs, 3);
hf = plot_kldiv_mat(med_kl_divs, chans, sprintf('Median over %d recordings', n_recs));
savefig(hf, fullfile(results_dir, 'res_figs', 'med_kl_div.fig'));

%%
function hf = plot_kldiv_mat(kl_div, chans, title_line2)

hf = figure;
sanePColor(kl_div);
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

ylabel('Channel i');
xlabel('Channel j');
end