% Do NMF and extract classes and transitions

%% Gather recordings from each day
sr_dirs = prepSR;
grid_rec_dir = fullfile(sr_dirs.results, 'grid_recs');

% These recs are the baseline ones that don't have BS
days = {
    '2019-03-07'
%     '2019-03-14'
%     '2019-04-25'
%     '2019-08-27'
    };

n_days = length(days);

input_s = struct('name', days, ...
                 'mt_res_in', cell(n_days, 1), ...
                 'nmf_res_out', cell(n_days, 1), ...
                 'xval_fig_dir', cell(n_days, 1));

for kD = 1:n_days
    % Build input to concat_and_nmf
    
    curr_day = days{kD};    
    input_s(kD).nmf_res_out = fullfile(grid_rec_dir, 'nmf_res', sprintf('nmf_res_%s.mat', curr_day));
    xval_dir = fullfile(grid_rec_dir, 'xval_figs', curr_day);
    
    if ~exist(xval_dir, 'dir') && ~mkdir(xval_dir)
        error('Could not create folder %s for cross-validation figures', xval_dir);
    end
    input_s(kD).xval_fig_dir = xval_dir;
    
    % find all mt_res files for this day
    res_entries = dir(fullfile(grid_rec_dir, 'mt_res', sprintf('mt_res_%s_*.mat', curr_day)));
    fns = sort({res_entries.name});
    input_s(kD).mt_res_in = fullfile(grid_rec_dir, 'mt_res', fns);
end

%% Do NMF
nmf_mfiles = concat_and_nmf(input_s);

%% Generate null model data
gen_null_model_data(nmf_mfiles);

%% Do KL divergence analysis
kl_div_subsets = struct('grid', 'G*', 'fork', 'F*');
kl_divergence_analysis(nmf_mfiles, kl_div_subsets);

%% Evaluate relationship b/w channel distance and KL divergence

kl_div_info = struct;
subsets = fieldnames(kl_div_subsets);
n_sets = length(subsets);
for kS = 1:n_sets
    for kl_field = ["distance", "kl_div", "kl_div_null"]
        kl_div_info.(subsets{kS}).(kl_field) = cell(n_days, 1);
    end
end

for kD = 1:n_days
    % Get channel locations from mt results
    mt_res_mfile = matfile(input_s(kD).mt_res_in{1});
    chan_locs = mt_res_mfile.chan_locs;   
    
    all_chan_names = mt_res_mfile.name;
    fork_chan_names = all_chan_names(strncmp(all_chan_names, 'F', 1));
    grid_chan_names = all_chan_names(strncmp(all_chan_names, 'G', 1));
    
    for kS = 1:n_sets
        set_name = subsets{kS};
        set_chans = nmf_mfiles{kD}.([set_name, '_chans']);
        
        in_grid = all(ismember(set_chans, grid_chan_names));
        in_fork = all(ismember(set_chans, fork_chan_names));
        
        if ~(in_grid || in_fork)
            warning('Cannot do distance analysis for channel set "%s" because its channels are not all from the same probe', set_name);
            kl_div_info = rmfield(kl_div_info, set_name);
            continue;
        end
        
        % identify the indices within the chan_locs struct that correspond to this channel subset
        if in_grid
            chan_loc_inds = cellfun(@(ch) find(strcmp(ch, grid_chan_names)), set_chans);            
            set_chan_locs = chan_locs.Probe1(chan_loc_inds, :); % Probe1 is the grid
        else
            chan_loc_inds = cellfun(@(ch) find(strcmp(ch, fork_chan_names)), set_chans);
            set_chan_locs = chan_locs.Probe2(chan_loc_inds, :); % Probe2 is the linear probe (fork)
        end
        
        dist_mat = squareform(pdist(set_chan_locs));
        kl_div_info.(set_name).distance{kD} = reshape(dist_mat, [], 1);
        
        % add in KL divergence (real and null)
        kl_div = nmf_mfiles{kD}.([set_name, '_divs']);
        kl_div_info.(set_name).kl_div{kD} = reshape(mean(kl_div, 3), [], 1);
        
        kl_div_null = nmf_mfiles{kD}.([set_name '_divs_null']);
        kl_div_info.(set_name).kl_div_null{kD} = reshape(mean(kl_div_null, 3), [], 1);
    end
end

save(fullfile(grid_rec_dir, 'kl_div_combined.mat'), '-struct', 'kl_div_info', '-v7.3');

%% Scatter-plot KL div and null KL div vs. distance, for each set. For now combine across all days.

% load in case we have to
if ~exist('kl_div_info', 'var')
    kl_div_info = load(fullfile(grid_rec_dir, 'kl_div_combined.mat'));
end

subsets = fieldnames(kl_div_info);
n_sets = length(subsets);
for kS = 1:n_sets
    distances = cell2mat(kl_div_info.(subsets{kS}).distance);
    kl_divs = cell2mat(kl_div_info.(subsets{kS}).kl_div);
    kl_divs_null = cell2mat(kl_div_info.(subsets{kS}).kl_div_null);
    
    hf = figure;
    scatter(distances, kl_divs, 'DisplayName', 'Real states');
    hold on;
    scatter(distances, kl_divs_null, '.', 'DisplayName', 'Null model');
    title(['KL divergences between NMF scores in ', subsets{kS}]);
    xlabel('Pair distance ({\mu}m)');
    ylabel('KL divergence (bits)');
    legend;
    
    savefig(hf, sprintf('%s_dist_DKL.fig', subsets{kS}));
end
