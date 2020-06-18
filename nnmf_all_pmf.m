% Cluster states in each recording (individually) using non-negative matrix factorization.

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
    
    % all channels should have the same nans
    pxx = res_mfile.pxx(1, 1);
    b_nan = any(isnan(pxx{1}));
    time_axis = time_axis(~b_nan);

    %% Preprocess data (manually run inner cell to redo)
    if ~isprop(res_mfile, 'pxx_rankord')
        %%
        % Here I am trying the rank-order "normalization" method.

        pp_options = struct;
        pp_options.in_name = 'pxx';
        pp_options.name = 'pxx_rankord';
        pp_options.freq_sm_type = 'med';
        pp_options.freq_sm_span = 10;
        pp_options.time_sm_type = 'exp';
        pp_options.time_sm_span = 120;
        pp_options.norm_type = 'rank';

        mt_preprocess(res_mfile, pp_options);
    
    end
    
    %% Get data for all channels and remove nans
    rec_data = res_mfile.pxx_rankord;
    rec_data = cellfun(@(d) d(:, ~b_nan), rec_data, 'uni', false);
    
    % allocate matrices to store U, V, and classifications
    U_all = cell(n_chans, 1);
    V_all = cell(n_chans, 1);
    U_pmf_all = cell(n_chans, 1);
    V_pmf_all = cell(n_chans, 1); % probability mass function, including "zero class"
    classes = zeros(length(time_axis), n_chans);
    
    
    %% Rescale and classify for each channel separately
    for kC = 1:n_chans
        
        % Test how many NMF components satisfactorily reduce error (manually run inner cell to redo)
        if ~isprop(res_mfile, ['nmf_comps_', chan_vnames{kC}])
            %%
            [n_comps, hfig] = nmf_ncomps_xval(rec_data{kC}(:, 1:20:end).', 1, 15, 0.01);

            figure(hfig);
            title(sprintf('NMF cross-validation (%s, %s)', rec, chans{kC}));        
            savefig(hfig, fullfile(results_dir, rec, ['nmf_xval_', chan_vnames{kC}, '.fig']));
            res_mfile.(['nmf_comps_', chan_vnames{kC}]) = n_comps;        
        end
        
        %% Do NMF on full dataset
        n_comps = res_mfile.(['nmf_comps_', chan_vnames{kC}]);
        [V, U, p] = sp_nnmf(rec_data{kC}.', n_comps, [], [], 500000);

        
        %% Sort comps according to peak frequency, and normalize power of U, then add all-0 component
        
        % sort
        [~, peak_freqinds] = max(U);
        % mean_freqs = freq_axis * (U ./ sum(U));
        [~, order] = sort(peak_freqinds);
        U = U(:, order);
        V = V(:, order);
        
        % normalize
        norm_factor = vecnorm(U);
        U = U ./ norm_factor;
        V = V .* norm_factor;
        
        % add all-0 component, normalize so sum of scores at each point is 1
        max_V_sum = max(sum(V, 2));
        V_pmf = V / max_V_sum;
        V_pmf = [V_pmf, 1 - sum(V_pmf, 2)];
        U_pmf = [U * max_V_sum, zeros(size(U, 1), 1)];

        U_pmf_all{kC} = U_pmf;
        V_pmf_all{kC} = V_pmf;
        
        %% plot components, component-space representation and reconstruction
        
        % plot components
        figure;
        sanePColor(1:n_comps, freq_axis, U, false, true);
        set(gca, 'YScale', 'log');
        xticks(1:n_comps);
        xlabel('Component #');
        ylabel('Frequency (Hz)');
        title(sprintf('NMF components (%s, %s)', chans{kC}, rec));

        % plot representation
        figure;
        sanePColor(time_axis, 1:n_comps, V.');
        xlabel('Time (s)');
        yticks(1:n_comps);
        ylabel('Component #');
        title(sprintf('NMF component representation of %s on %s', chans{kC}, rec));
        
        % plot reconstruction
        figure;
        sanePColor(time_axis, freq_axis, U * V', false, true);
        set(gca, 'YScale', 'log');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(sprintf('NMF reconstruction of %s on %s (rank)', chans{kC}, rec));

    end  
    %% Save to the matfile
    res_mfile.pmf_U = U_pmf_all;
    res_mfile.pmf_V = V_pmf_all;
end

function [n_comps, hfig] = nmf_ncomps_xval(data, min_comps, max_comps, pow_thresh)
% Do cross-validation to test how many NMF components to use.
% Inputs:
%   data:        the non-negative data matrix to decompose
%   min_comps:   smallest # of components to try
%   max_comps:   largest # of components to try
%   pow_thresh:  fraction of explained input power necessary to justify adding each component
% Outputs:
%   n_comps:     best # of components
%   hfig:        handle to figure showing cross-validation results

n_comp_options = min_comps:max_comps;
BCV_err = nnmf_k_ver(data, round(size(data) / 5), 5, min_comps, max_comps, [], 500000);

total_power = norm(data, 'fro');
explained = 1 - BCV_err / total_power;

hfig = figure;
scatter(reshape(repmat(n_comp_options, size(explained, 1), 1), [], 1), explained(:), 'k', 'filled');
hold on;
plot(n_comp_options, mean(explained), 'b', 'LineWidth', 1);
xticks(n_comp_options);
xlabel('Number of NMF components');
ylabel('Fraction of explained power');
title('NMF cross-validation');

% use all components that explain at least 1% of power
n_comps = find(diff(mean(explained)) < pow_thresh, 1);

if isempty(n_comps)
    n_comps = n_comp_options(end);
    warning('All components explain > 1% of power - consider trying more.');
end

ylims = get(gca, 'YLim');
h = plot([1 1] * (n_comps + 0.5), ylims, 'r--');
ylim(ylims);
legend(h, '1% power threshold', 'Location', 'southeast');

end