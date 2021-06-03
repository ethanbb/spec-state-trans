function vals_bytype = categorize_pair_data(mats_over_days, chan_names)
% Take some matrices of data over days (in 3rd dimension) and return struct w/ combined
% data vectors, with nans removed, in these mutually exclusive categories:
% 'SameChannel', 'SameRegionL4', 'SameRegionNonL4', 'CrossRegionL4', 'CrossRegionNonL4'
% Any pair of a channel with itself is 'SameChannel'
% Any pair where either channel is labeled as L4 is in the corresponding L4 category
% Other pairs are in the corresponding NonL4 category depending on whether the region is the same.
% chan_names should be a cell array of underscore-separated channel names (e.g. 'V1R_Sup1')
% which is the same length as the first two dimensions of mats_over_days.

chan_names = string(chan_names(:));
n_chans = length(chan_names);
chan_name_parts = split(chan_names, "_");
regions = chan_name_parts(:, 1);
layers = chan_name_parts(:, 2);

layer_isl4 = layers == "L4";
pair_hasl4 = layer_isl4 | layer_isl4';
same_region = regions == regions';
same_channel = logical(eye(n_chans));

% to broadcast along later dimensions
if ismatrix(mats_over_days)
    get_masked_vals = @(mask) mats_over_days(mask);
else
    broadcast_dims = size(mats_over_days, 3:ndims(mats_over_days));
    get_masked_vals = @(mask) mats_over_days(repmat(mask, [1, 1, broadcast_dims]));   
end

vals_bytype.SameChannel = get_masked_vals(same_channel);
vals_bytype.SameRegionL4 = get_masked_vals(~same_channel & pair_hasl4 & same_region);
vals_bytype.SameRegionNonL4 = get_masked_vals(~same_channel & ~pair_hasl4 & same_region);
vals_bytype.CrossRegionL4 = get_masked_vals(~same_channel & pair_hasl4 & ~same_region);
vals_bytype.CrossRegionNonL4 = get_masked_vals(~same_channel & ~pair_hasl4 & ~same_region);
vals_bytype.SameRegion = [vals_bytype.SameRegionL4; vals_bytype.SameRegionNonL4];
vals_bytype.CrossRegion = [vals_bytype.CrossRegionL4; vals_bytype.CrossRegionNonL4];

vals_bytype = structfun(@(v) v(~isnan(v)), vals_bytype, 'uni', false);

end
