sr_dirs = prepSR;

[exp_info, info_s_all] = gather_exp_info;
run_dates = {input_s_all.name}';
m1v1_dates = input_s_all(1).days;
bilatv1_dates = input_s_all(2).days;

nmf_res_files = fullfile(run_dates, 'nmf_res.mat');

[characteristic_sortinds, top_100_inds] = cellfun(@(f) ...
    sort_class_data_by_typicality(f, [], false, false), nmf_res_files, 'uni', false);

save('characteristic_inds.mat', 'run_dates', 'characteristic_sortinds', 'top_100_inds', ...
    'm1v1_dates', 'bilatv1_dates');
