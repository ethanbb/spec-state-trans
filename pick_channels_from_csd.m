% For all recording days, use CSD analysis and diagram of rat cortical layers
% to pick channels to analyze. Aiming for 1 superficial (layer 2/3), 1 layer 4, and 1 layer 5.

% Get recording days

sr_dirs = prepSR;

days = {
    '2020-01-30'
    '2020-01-31'
    '2020-02-06'
    '2020-03-05'
    '2020-03-06'
    '2020-03-10'
    '2020-03-11'
    '2020-10-26'
    '2020-10-27'
    '2020-10-28'
    '2020-10-29'
    };
n_days = length(days);

csd_dirs = fullfile(sr_dirs.results, days);

for kD = 1:n_days
    uiwait(pick_csd_channels(csd_dirs{kD}));
end


