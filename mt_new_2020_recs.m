% Do CSD and multitaper on the October 2020 recordings

sr_dirs = prepSR;

rec_names = {
    '2020-10-26_11-30-00'
    '2020-10-26_13-31-00'
    '2020-10-26_15-33-00'
    '2020-10-27_13-04-00'
    '2020-10-27_15-05-00'
    '2020-10-28_12-29-00'
    '2020-10-28_14-31-00'
    '2020-10-28_16-45-00'
    '2020-10-29_12-05-00'
    '2020-10-29_14-06-00'
    '2020-10-29_16-07-00'
    };

temp = split(rec_names, '_');
rec_dates = unique(temp(:, 1));

%% Do CSD

probe_s = struct('Probe1', 'M1', 'Probe2', 'V1');
for kD = 1:length(rec_dates)
    plot_csd(rec_dates{kD}, probe_s);
end
