function nmf_mfiles = concat_and_nmf(input_s, pp_options_s)
% Concatenate datasets and do non-negative matrix factorization
% to reduce dimensionality.
%
% Syntax: concat_and_nmf(input_s, savedir)
% Inputs:
%   input_s: struct array, each element corresponds to an output file and has keys:
%       name: a user-friendly identifier for this set of inputs (e.g. the day)
%       mt_res_in: cell array of mt_res file paths
%       nmf_res_out: output file path
%       xval_fig_dir: place to save NMF cross-validation figures
%
%   pp_options_s: optional struct of options for mt_preprocess to override defaults. See comments on 'opts'
%                 in mt_preprocess for details. Defaults are in the "Preprocess data" cell below.

nmf_mfiles = arrayfun(@(input) matfile(input.nmf_res_out, 'Writable', true), input_s, 'uni', false);

for kO = 1:length(input_s)
    input_info = input_s(kO);
    curr_name = input_info.name;
    save_mfile = nmf_mfiles{kO};

    res_mfiles = cellfun(@(p) matfile(p, 'Writable', true), ...
        input_info.mt_res_in, 'uni', false);
    n_files = length(res_mfiles);
    
    chans = res_mfiles{1}.name;
    chan_vnames = cellfun(@matlab.lang.makeValidName, chans, 'uni', false); % valid for variables etc.
    n_chans = length(chans);

	mt_opts = res_mfiles{1}.options;
    freq_axis = res_mfiles{1}.freq_grid;
    
    % preprocessing defaults
    pp_options = struct;
    pp_options.in_name = 'pxx';      % input field name
    pp_options.name = 'pxx_rankord'; % output field name
    pp_options.freq_sm_type = 'med'; % median smooth over frequencies for
    pp_options.freq_sm_span = 10;    % 10 frequency bins
    pp_options.time_sm_type = 'exp'; % exponential smooth in time for
    pp_options.time_sm_span = 120;   % 120 seconds
    pp_options.norm_type = 'rank';   % Rank-order "normalization"
    
    % override by input if provided
    if exist('pp_options_s', 'var') && ~isempty(pp_options_s)
        pp_fields = fieldnames(pp_options_s);
        for kF = 1:length(pp_fields)
            pp_options.(pp_fields{kF}) = pp_options_s.(pp_fields{kF});
        end
    end
    
    %% Preprocess data (manually run inner cell to redo)
    for kF = 1:n_files
        if ~isprop(res_mfiles{kF}, pp_options.name)
            %%
            mt_preprocess(res_mfiles{kF}, pp_options);
        end
    end
    %% Get data for all channels and concatenate good segments
    
    % Also create a time axis for all recordings combined, with the second one picking
    % up where the first one ends etc.
    
    time_axes = cell(1, n_files);
    rec_data = repmat({cell(1, n_files)}, n_chans, 1);
    
    last_endpoint = 0; % End of time covered by previous recording
    win_length = mt_opts.window; % Add half of this after last timepoint to account for window length
    for kF = 1:n_files
        pxx_processed = res_mfiles{kF}.(pp_options.name);
        seg_windows = res_mfiles{kF}.seg_windows;
        good_windows_cat = cell2mat(seg_windows);

        for kC = 1:n_chans
            rec_data{kC}{kF} = pxx_processed{kC}(:, good_windows_cat);
        end
        
        time_axes{kF} = res_mfiles{kF}.time_grid + last_endpoint;
        last_endpoint = time_axes{kF}(end) + win_length/2;
        time_axes{kF} = time_axes{kF}(good_windows_cat);
    end
    time_axis = cell2mat(time_axes);
    rec_data = cellfun(@cell2mat, rec_data, 'uni', false); % concatenate each channel individually

    % allocate matrices to store U, V, and classifications
    % do it twice in order to be able to compare each channel with itself.
    n_reps = 2;
    U_all = cell(n_reps, 1);
    V_all = cell(n_reps, 1);
    classes_all = cell(n_reps, 1);
    trans_all = cell(n_reps, 1);

    if ~isprop(save_mfile, 'nmf_comps')
        need_ncomps = true(n_chans, 1);
        save_mfile.nmf_comps = nan(n_chans, 1);
    else
        need_ncomps = isnan(save_mfile.nmf_comps);
    end
    nmf_comps = save_mfile.nmf_comps;
      
%% Rescale and classify for each channel separately

    for kR = 1:n_reps
        U_rep = cell(n_chans, 1);
        V_rep = cell(n_chans, 1);
        classes_rep = cell(n_chans, 1);
        trans_rep = cell(n_chans, 1);
        
        xval_fig_dir = input_info.xval_fig_dir;
        
        parfor kC = 1:n_chans

            % Test how many NMF components satisfactorily reduce error (manually run inner cell to redo)
            if need_ncomps(kC)
                %%
                [n_comps, hfig] = nmf_ncomps_xval(rec_data{kC}(:, 1:20:end).', 1, 15, 0.01);

                figure(hfig);
                title(sprintf('NMF cross-validation (%s, %s)', curr_name, chans{kC}));        
                savefig(hfig, fullfile(xval_fig_dir, ['nmf_xval_', chan_vnames{kC}, '.fig']));
                nmf_comps(kC) = n_comps;
                need_ncomps(kC) = false;

                %%
                close(hfig);
            end
        
            %% Do NMF on full dataset
            n_comps = nmf_comps(kC);
            [V, U] = sp_nnmf(rec_data{kC}.', n_comps, [], [], 500000);
            % Version that tries to make scores sparser, but it's not necessary
            %[U, V] = kim_park_snmf(rec_data{kC}, n_comps, 0, 0.4, [], [], 500000);

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
            % copied from pmtk3's countTransitions
            trans = accumarray([classes(1:end-1), classes(2:end)], 1, [n_comps, n_comps]);

            U_rep{kC} = U;
            V_rep{kC} = V;
            classes_rep{kC} = classes;
            trans_rep{kC} = trans;
                
            %% plot components, component-space representation and reconstruction

%             % plot components
%             figure;
%             sanePColor(1:n_comps, freq_axis, U, false, true);
%             set(gca, 'YScale', 'log');
%             xticks(1:n_comps);
%             xlabel('Component #');
%             ylabel('Frequency (Hz)');
%             title(sprintf('NMF components (%s, %s)', chans{kC}, curr_name));
% 
%             % plot representation
%             figure;
%             sanePColor(time_axis, 1:n_comps, V.');
%             xlabel('Time (s)');
%             yticks(1:n_comps);
%             ylabel('Component #');
%             title(sprintf('NMF component representation of %s on %s', chans{kC}, curr_name));
%             
%             % plot reconstruction
%             figure;
%             sanePColor(time_axis, freq_axis, U * V', false, true);
%             set(gca, 'YScale', 'log');
%             xlabel('Time (s)');
%             ylabel('Frequency (Hz)');
%             title(sprintf('NMF reconstruction of %s on %s (rank)', chans{kC}, curr_name));

        end
        
        U_all{kR} = U_rep;
        V_all{kR} = V_rep;
        classes_all{kR} = classes_rep;
        trans_all{kR} = trans_rep;
    end  
    %% Save to the matfile
    save_mfile.run_name = curr_name;
    save_mfile.nmf_comps = nmf_comps;
    save_mfile.chan_names = chans;
    save_mfile.nmf_U = U_all;
    save_mfile.nmf_V = V_all;
    save_mfile.nmf_classes = classes_all;
    save_mfile.nmf_transitions = trans_all;
    save_mfile.time_axis = time_axis;
    save_mfile.freq_axis = freq_axis;
    
end

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

hfig = figure('Visible', 'off');
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
