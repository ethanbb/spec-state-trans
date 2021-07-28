function vals_bytype = categorize_pair_data(mats_over_days, chan_names)
% Take some matrices of data over days (in 3rd dimension) and return struct w/ combined
% data vectors, with nans removed, in these mutually exclusive categories:
% 'SameChannel', 'SameRegionL4', 'SameRegionNonL4', 'CrossRegionL4', 'CrossRegionNonL4'
% Any pair of a channel with itself is 'SameChannel'
% Any pair where either channel is labeled as L4 is in the corresponding L4 category
% Other pairs are in the corresponding NonL4 category depending on whether the region is the same.
% chan_names should be a cell array of underscore-separated channel names (e.g. 'V1R_Sup1')
% which is the same length as the first two dimensions of mats_over_days.
%
% This function has been mostly outsourced to get_symmetric_matrix_masks
% but kept for old code.

s_masks = util.get_symmetric_matrix_masks(chan_names);

% to broadcast along later dimensions
if ismatrix(mats_over_days)
    get_masked_vals = @(mask) mats_over_days(mask);
else
    broadcast_dims = size(mats_over_days, 3:ndims(mats_over_days));
    get_masked_vals = @(mask) mats_over_days(repmat(mask, [1, 1, broadcast_dims]));   
end

vals_bytype = structfun(@(mask) get_masked_vals(mask), s_masks, 'uni', false);
vals_bytype = structfun(@(v) v(~isnan(v)), vals_bytype, 'uni', false);

end
