% Fix snippit files to make sure they are all aligned to the correct stimulus onsets

recs_to_fix = {
    '2021-01-29_10-39-00'
    '2021-01-29_12-09-00'
    '2021-01-29_16-53-00'
    '2021-01-31_09-52-00'
    '2021-01-31_11-21-00'
    '2021-01-31_12-39-00'
    '2021-02-01_12-10-00'
    '2021-02-01_13-48-00'
    '2021-02-02_10-42-00'
    '2021-02-02_12-13-00'
    '2021-02-02_13-43-00'
    '2021-02-02_15-14-00'
    };

sr_dirs = prepSR;
p = addpath('Z:\code\Brenna_code\preprocessing');

% parameters for extractSnippits
before = 1;
l = 3;
finalSampR = 1000;
startBaseline = 0;
lfp_field = 'meanSubFullTrace';

for kR = 1:length(recs_to_fix)
    rec = recs_to_fix{kR};
    lfp_data_s = load(fullfile(sr_dirs.processed_lfp, sprintf('meanSub_%s.mat', rec)), lfp_field);
    lfp_data = lfp_data_s.(lfp_field);
    snips_mfile = matfile(fullfile(sr_dirs.snippits, sprintf('snips_%s.mat', rec)), 'Writable', true);
    
    starts = snips_mfile.allStartTimes;
    new_snips1 = extractSnippets_Plexon_BPS(lfp_data(1:64, :), starts{2}, before, l, finalSampR, startBaseline);
    new_snips2 = extractSnippets_Plexon_BPS(lfp_data(65:end, :), starts{1}, before, l, finalSampR, startBaseline);
    
    n = min(size(new_snips1, 2), size(new_snips2, 2));
    new_snips = vertcat(new_snips1(:, 1:n, :), new_snips2(:, 1:n, :));
    snips_mfile.dataSnippits = new_snips;
end

path(p);
