% Cluster states for each channel on each day using non-negative matrix factorization.

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

input_s = struct('name', days, ...
                 'mt_res_in', cell(n_days, 1), ...
                 'nmf_res_out', cell(n_days, 1), ...
                 'xval_fig_dir', cell(n_days, 1));

for kD = 1:n_days
    % Build input to concat_and_nmf
    
    curr_day = days{kD};
    input_s(kD).nmf_res_out = fullfile(results_dir, curr_day, 'nmf_res.mat');
    input_s(kD).xval_fig_dir = fullfile(results_dir, curr_day);
    
    % find all "layers" results files under this day
    res_fn = 'mt_res_layers.mat';
    res_entries = dir(fullfile(results_dir, curr_day, '*', res_fn));
    dirs = sort({res_entries.folder});
    input_s(kD).mt_res_in = fullfile(dirs, res_fn);
end

concat_and_nmf(input_s);
