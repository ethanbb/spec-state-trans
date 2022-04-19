sr_dirs = prepSR;

run_dates = {
    '2020-10-26'
    '2020-10-28'
    '2020-10-29'
    '2021-01-27'
    '2021-01-29'
    '2021-01-31'
    '2021-02-02'
    };

m1v1_dates = run_dates(1:3);
bilatv1_dates = run_dates(4:7);

nmf_res_files = fullfile(run_dates, 'nmf_res.mat');

[characteristic_sortinds, top_100_inds] = cellfun(@(f) ...
    sort_class_data_by_typicality(f, [], false, false), nmf_res_files, 'uni', false);

save('characteristic_inds.mat', 'run_dates', 'characteristic_sortinds', 'top_100_inds', ...
    'm1v1_dates', 'bilatv1_dates');
