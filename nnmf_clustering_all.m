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
    rec = recs{kR};
    
    %% Get dataset info
    res_path = fullfile(results_dir, rec, 'mt_res.mat');
    res_mfile = matfile(res_path, 'Writable', true);
    
    mt_opts = res_mfile.options; % necessary due to MatFile quirk
    Fw = 1 / mt_opts.winstep; % "window rate"
    
    chans = res_mfile.name;
    chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false);
    n_chans = numel(chans);
    assert(n_chans == 2, 'Deal with # of channels in %s', rec);
    
    %% Preprocess data (skippable if already done)
    % Here I am trying the rank-order "normalization" method.
    
    pp_options = struct;
    pp_options.in_name = 'pxx';
    pp_options.name = 'pxx_rankord';
    pp_options.freq_sm_type = 'med';
    pp_options.freq_sm_span = 10;
    pp_options.time_sm_type = 'exp';
    pp_options.time_sm_span = 60;
    pp_options.norm_type = 'rank';
    
    mt_preprocess(res_mfile, pp_options);
    
    %% Get data for all channels
    rec_data = res_mfile.pxx_rankord;
    freq_axis = res_mfile.freq_grid;
    time_axis = res_mfile.time_grid;
    
    %% Process each channel separately
    for kC = 1:n_chans
        %% Remove NaNs
        chan_data = rec_data{kC};
        b_nan = any(isnan(chan_data));
        chan_data_nonan = chan_data(:, ~b_nan);
        
        %% Test how many NMF components satisfactorily reduce error (skippable)
        n_comp_options = 1:15;
        probe_data = chan_data_nonan(:, 1:20:end).';
        BCV_err = nnmf_k_ver(probe_data, round(size(probe_data) / 10), 10, ...
            n_comp_options(1), n_comp_options(end), [], 500000);
        
        total_power = norm(probe_data, 'fro');
        explained = 1 - BCV_err / total_power;
        
        figure;
        scatter(reshape(repmat(n_comp_options, size(explained, 1), 1), [], 1), ...
            explained(:), 'k', 'filled');
        hold on;
        plot(n_comp_options, mean(explained), 'b', 'LineWidth', 1);
        xticks(n_comp_options);
        xlabel('Number of NMF components');
        ylabel('Fraction of explained power');
        title(sprintf('NMF cross-validation (%s, %s)', rec, chans{kC}));
        
        % use all components that explain at least 1% of power
        n_comps = find(diff(mean(explained)) < 0.01, 1);
        if isempty(n_comps)
            n_comps = n_comp_options(end);
            warning('All components explain > 1% of power - consider trying more.');
        end
        ylims = get(gca, 'YLim');
        h = plot([1 1] * (n_comps + 0.5), ylims, 'r--');
        legend(h, '1% power threshold', 'Location', 'southeast');
        
        savefig(fullfile(results_dir, rec, ['nmf_xval_', chan_vnames{kC}, '.fig']));
        res_mfile.(['nmf_comps_', chan_vnames{kC}]) = n_comps;
        
        %% Do NMF on full dataset
        n_comps = res_mfile.(['nmf_comps_', chan_vnames{kC}]);
        [U, V, p] = sp_nnmf(chan_data_nonan.', n_comps, [], [], 500000);
        
        %% Sort comps according to peak frequency and plot
        [~, colmax] = max(V);
        [~, order] = sort(colmax);
        U = U(:, order);
        V = V(:, order);
        comp_axis = 1:n_comps;
        
        figure;
        sanePColor(comp_axis, freq_axis, V, false, true);
        set(gca, 'YScale', 'log');
        xticks(1:n_comps);
        xlabel('Component #');
        ylabel('Frequency (Hz)');
        title('NMF components');

        %% Visualize NMF result along time domain
        U_full = nan(length(time_axis), n_comps);
        U_full(~b_nan, :) = U;
        
        % plot component-space representation
        figure;
        sanePColor(time_axis, comp_axis, U_full.');
        xlabel('Time (s)');
        ylabel('Component #');
        title(sprintf('NMF component representation of %s on %s', chans{kC}, rec));
        
        % plot reconstruction
        figure;
        sanePColor(time_axis, freq_axis, (U_full * V')', false, true);
        set(gca, 'YScale', 'log');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(sprintf('NMF reconstruction of %s on %s (rank)', chans{kC}, rec));
    end
end
