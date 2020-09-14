% Do NMF and extract classes and transitions

prepSR;

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
    input_s(kD).nmf_res_out = fullfile(results_dir, 'grid_recs', 'nmf_res', ...
        sprintf('nmf_res_%s.mat', curr_day));
    input_s(kD).xval_fig_dir = fullfile(results_dir, 'grid_recs', 'xval_figs', curr_day);
    
    if ~mkdir(input_s(kD).xval_fig_dir)
        error('Could not create folder %s for cross-validation figures', input_s(kD).xval_fig_dir);
    end
    
    % find all mt_res files for this day
    res_entries = dir(fullfile(results_dir, 'grid_recs', 'mt_res', sprintf('mt_res_%s_*.mat', curr_day)));
    fns = sort({res_entries.name});
    input_s(kD).mt_res_in = fullfile(results_dir, 'grid_recs', 'mt_res', fns);
end

concat_and_nmf(input_s);
