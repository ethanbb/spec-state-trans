function bs_segments = find_likely_bs(data_mfile, chans_to_look_at, bs_dur_thresh)
% Given a matfile of LFP data, find segments that might have burst suppression
% (to be examined manually). The test is pretty sensitive; if this doesn't find
% anything, it's probably safe to assume that there's no burst suppression.
% chans_to_look_at is an array of channel indices (after ordering) to process and merge segments of;
% if omitted, defaults to the first and last channel (to get one for each probe in a 2-probe setup).
% bs_dur_thresh is the min duration at near isoelectric potential to consider burst suppression
% initially (before defragmenting segments). Defaults to 2 (seconds).

if ischar(data_mfile)
    data_mfile = matfile(data_mfile);
end

if ~exist('chans_to_look_at', 'var') || isempty(chans_to_look_at)
    n = size(data_mfile, 'meanSubFullTrace', 1);
    chans_to_look_at = [1, n];
end

if ~exist('bs_dur_thresh', 'var') || isempty(bs_dur_thresh)
    bs_dur_thresh = 2;
end

Fs = data_mfile.finalSampR;
lfp = organize_lfp(data_mfile, chans_to_look_at);
b_bs_mat = find_bs(lfp', Fs, bs_dur_thresh);
b_bs = any(b_bs_mat, 2);
b_bs_combined = combine_bs_segments(b_bs, Fs);

% convert to start and end points in seconds
bs_segments = find_segments(b_bs_combined, true, Fs);
% round 
bs_segments(:, 1) = floor(bs_segments(:, 1) / Fs);
bs_segments(:, 2) = min(ceil(bs_segments(:, 2) / Fs), length(b_bs_combined) / Fs);

end

function b_bs = find_bs(lfp_mat, Fs, bs_dur_thresh)
% finds burst suppression periods in each column of lfp_mat

% make an object to execute moving root-mean-square (BS should have low RMS)
rms_obj = dsp.MovingRMS('Method', 'Exponential weighting', 'ForgettingFactor', 0.99^(1000/Fs));
moving_rms = rms_obj(lfp_mat);

% find "suppression" periods with findpeaks
b_bs = false(size(moving_rms));
for kC = 1:size(moving_rms, 2)
    [n, binedges] = histcounts(moving_rms(:, kC));
    bincenters = mean([binedges(1:end-1); binedges(2:end)]);
    [~, locs, ~, proms] = findpeaks(-n, bincenters);
    [~, imaxprom] = max(proms);
    thresh = locs(imaxprom);
    
    b_bs = moving_rms(:, kC) < thresh;
end

% find non-suppression periods by thresholding RMS and smoothing. points that remain 0 are suppression.
% b_bs = double(moving_rms < 0.1);
sm_b_bs = smoothdata(b_bs, 'movmean', round(bs_dur_thresh * Fs));
b_bs = sm_b_bs == 1;
end

function b_bs = combine_bs_segments(b_bs, Fs)
% mark short blocks of good data between suppression blocks as also suppression
[good_segs, good_durs] = find_segments(b_bs, false, Fs);
mark_supp_inds = find(good_durs < 60);

for ind = mark_supp_inds(:)'
    b_bs(good_segs(ind, 1):good_segs(ind, 2)) = true;
end

% now filter out bs segments of < 20 seconds that are isolated
[bs_segs, bs_durs] = find_segments(b_bs, true, Fs);
mark_good_inds = find(bs_durs < 20);

for ind = mark_good_inds(:)'
    b_bs(bs_segs(ind, 1):bs_segs(ind, 2)) = false;
end

% one more round, removing good segments we consider too short
[good_segs, good_durs] = find_segments(b_bs, false, Fs);
mark_supp_inds = find(good_durs < 60 * 5);

for ind = mark_supp_inds(:)'
    b_bs(good_segs(ind, 1):good_segs(ind, 2)) = true;
end

end

function [segs, durs] = find_segments(b_bs, find_bs, Fs)
% Find segments of either burst suppression or non-burst suppression (depending on find_bs) in b_bs.
% Returns an n_segs x 2 matrix of start/end samples.
% durs is a vector of the duration of each segment in seconds.

b_bs_aug = [~find_bs; b_bs(:); ~find_bs];
if find_bs
    starts = find(diff(b_bs_aug) > 0);
    ends = find(diff(b_bs_aug) < 0) - 1;
else
    starts = find(diff(b_bs_aug) < 0);
    ends = find(diff(b_bs_aug) > 0) - 1;
end
segs = [starts, ends];
durs = diff(segs, [], 2) / Fs;

end
