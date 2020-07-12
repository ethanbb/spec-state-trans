% Cluster states for each channel on each day using non-negative matrix factorization.

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

%% Loop through days
for kD = 1:n_days
    %% Get dataset info
    
    day = days{kD};
    save_mfile = matfile(fullfile(results_dir, day, 'nmf_res.mat'), 'Writable', true);
    
    % find all "layers" results files under this day
    res_fn = 'mt_res_layers.mat';
    res_entries = dir(fullfile(results_dir, day, '*', res_fn));
    dirs = sort({res_entries.folder});
    res_paths = fullfile(dirs, res_fn);
    res_mfiles = cellfun(@(p) matfile(p, 'Writable', true), res_paths, 'uni', false);
    n_files = length(res_mfiles);
    
    chans = res_mfiles{1}.name;
    chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false); % valid for variables etc.
    n_chans = length(chans);
    
    % make sure all recordings use the same channels
    mt_opts = res_mfiles{1}.options;
    chan_inds = mt_opts.chans;
    
    for kF = 2:n_files
        this_opts = res_mfiles{kF}.options;
        this_chan_inds = this_opts.chans;
        assert(length(this_chan_inds) == length(chan_inds) && all(this_chan_inds == chan_inds), ...
            'Channels do not match between recordings on %s', day);
    end

    freq_axis = res_mfiles{1}.freq_grid;

    %% Preprocess data (manually run inner cell to redo)
    for kF = 1:n_files
        if ~isprop(res_mfiles{kF}, 'pxx_rankord')
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

            mt_preprocess(res_mfiles{kF}, pp_options);

        end
    end
    %% Get data for all channels and remove nans
    
    % Also create a time axis for all recordings combined, with the second one picking
    % up where the first one ends etc.
    
    time_axes = cell(1, n_files);
    b_nans = cell(1, n_files);
    rec_data = repmat({cell(1, n_files)}, n_chans, 1);
    
    last_endpoint = 0; % End of time covered by previous recording
    win_length = mt_opts.window; % Add half of this after last timepoint to account for window length
    for kF = 1:n_files
        pxx_rankord = res_mfiles{kF}.pxx_rankord;
        b_nans{kF} = all(isnan(pxx_rankord{1}));

        for kC = 1:n_chans
            rec_data{kC}{kF} = pxx_rankord{kC}(:, ~b_nans{kF});
        end
        
        time_axes{kF} = res_mfiles{kF}.time_grid + last_endpoint;
        last_endpoint = time_axes{kF}(end) + win_length/2;
        time_axes{kF} = time_axes{kF}(~b_nans{kF});
    end
    time_axis = cell2mat(time_axes);
    rec_data = cellfun(@cell2mat, rec_data, 'uni', false); % concatenate each channel individually

    % allocate matrices to store U, V, and classifications
    % do it twice in order to be able to compare each channel with itself.
    n_reps = 2;
    U_all = repmat({cell(n_chans, 1)}, n_reps, 1);
    V_all = repmat({cell(n_chans, 1)}, n_reps, 1);
    classes_all = repmat({cell(n_chans, 1)}, n_reps, 1);
    trans_all = repmat({cell(n_chans, 1)}, n_reps, 1);

    if ~isprop(save_mfile, 'nmf_comps')
        need_ncomps = true(n_chans, 1);
        save_mfile.nmf_comps = nan(n_chans, 1);
    else
        need_ncomps = isnan(save_mfile.nmf_comps);
    end
      
%% Rescale and classify for each channel separately

    for kC = 1:n_chans
        
        % Test how many NMF components satisfactorily reduce error (manually run inner cell to redo)
        if need_ncomps(kC)
            %%
            [n_comps, hfig] = nmf_ncomps_xval(rec_data{kC}(:, 1:20:end).', 1, 15, 0.01);

            figure(hfig);
            title(sprintf('NMF cross-validation (%s, %s)', day, chans{kC}));        
            savefig(hfig, fullfile(results_dir, day, ['nmf_xval_', chan_vnames{kC}, '.fig']));
            save_mfile.nmf_comps(kC, 1) = n_comps;
            
            %%
            close(hfig);
        end
        
        %% Do NMF on full dataset
        for kR = 1:n_reps
            n_comps = save_mfile.nmf_comps(kC, 1);
            [V, U, p] = sp_nnmf(rec_data{kC}.', n_comps, [], [], 500000);

            %% Sort comps according to peak frequency, and normalize power of U

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
            
            % get most likely "class"
            [~, classes] = max(V, [], 2);
            
            % count transitions between classes
            trans = countTransitions(classes, n_comps);
 
            U_all{kR}{kC} = U;
            V_all{kR}{kC} = V;
            classes_all{kR}{kC} = classes;
            trans_all{kR}{kC} = trans;
                
            %% plot components, component-space representation and reconstruction

%             % plot components
%             figure;
%             sanePColor(1:n_comps, freq_axis, U, false, true);
%             set(gca, 'YScale', 'log');
%             xticks(1:n_comps);
%             xlabel('Component #');
%             ylabel('Frequency (Hz)');
%             title(sprintf('NMF components (%s, %s)', chans{kC}, rec));
% 
%             % plot representation
%             figure;
%             sanePColor(time_axis, 1:n_comps, V.');
%             xlabel('Time (s)');
%             yticks(1:n_comps);
%             ylabel('Component #');
%             title(sprintf('NMF component representation of %s on %s', chans{kC}, rec));
%             
%             % plot reconstruction
%             figure;
%             sanePColor(time_axis, freq_axis, U * V', false, true);
%             set(gca, 'YScale', 'log');
%             xlabel('Time (s)');
%             ylabel('Frequency (Hz)');
%             title(sprintf('NMF reconstruction of %s on %s (rank)', chans{kC}, rec));

        end
    end  
    %% Save to the matfile
    save_mfile.nmf_U = U_all;
    save_mfile.nmf_V = V_all;
    save_mfile.nmf_classes = classes_all;
    save_mfile.nmf_transitions = trans_all;
    save_mfile.time_axis = time_axis;
    save_mfile.freq_axis = freq_axis;
    
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