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
% kl_divergence_analysis(nmf_mfiles);
