function [csd, valid_chans] = calc_kernel_csd(chan_data, bad_chans, b_output_interp, b_trim)
% Compute a 2nd spatial derivative CSD over the channels of chan_data.
% bad_chans is a vector of channel indices that should not be included.
% b_output_interp: values for bad channels are interpolated before doing the CSD. If this flag is
%                  false (default), replace the values on the bad channels with nans in the output.
% b_trim: (false) whether to remove first and last channels that are not at least half a kernel's
%         width from the edge. If false, these channels come back as nans.
%
% Output valid_chans is just those that are covered by the kernel (without regard for bad chans)

% defaults
b_output_interp = exist('b_output_interp', 'var') && b_output_interp;
b_trim = exist('b_trim', 'var') && b_trim;

% prep constants
num_chans = size(chan_data, 1);
sigma = 2;
support = -ceil(3*sigma):ceil(3*sigma);
n_edge = (length(support)-1) / 2;
valid_chans = 1+n_edge:num_chans-n_edge;

% hadle pathological case separately since it breaks interpolation
good_chans = setdiff(1:num_chans, bad_chans);
if isempty(good_chans)
    csd = nan(size(chan_data));
    if b_trim
        csd = csd(valid_chans, :);
    end
    return;
end

kernel = (-(sigma^2 - support.^2)/(sigma^5*sqrt(2*pi))) .* ...
    exp(-support.^2 / (2*sigma^2));

% interpolate over bad channels
chan_data(bad_chans, :) = interp1(good_chans, chan_data(good_chans, :), bad_chans);
csd = convn(chan_data, kernel.', 'valid');

if ~b_output_interp
    valid_bad_chans = intersect(bad_chans - n_edge, valid_chans);
    csd(valid_bad_chans, :) = nan;
end

if ~b_trim
    csd_valid = csd;
    csd = nan(size(chan_data));
    csd(valid_chans, :) = csd_valid;
end

end
