function s_masks = get_symmetric_matrix_masks(chan_names)
% Given the channel names for the rows and columns of a square matrix of pairwise statistics,
% make a struct of binary masks to select various sets of pairs.
% The matrix is assumed to be symmetric, and with the exception of SameChannel (which is
% the diagonal), all masks take values only from the lower triangle of the matrix.

chan_names = string(chan_names(:));
n_chans = length(chan_names);
lower_triangle = 1:n_chans < (1:n_chans)';

chan_name_parts = split(chan_names, "_");
regions = chan_name_parts(:, 1);
layers = chan_name_parts(:, 2);

layer_isl4 = layers == "L4";
pair_hasl4 = layer_isl4 | layer_isl4';
same_region = regions == regions';

s_masks.SameChannel = logical(eye(n_chans));
s_masks.SameRegion = lower_triangle & same_region;
s_masks.CrossRegion = lower_triangle & ~same_region;

s_masks.SameRegionL4 = lower_triangle & pair_hasl4 & same_region;
s_masks.SameRegionNonL4 = lower_triangle & ~pair_hasl4 & same_region;
s_masks.CrossRegionL4 = lower_triangle & pair_hasl4 & ~same_region;
s_masks.CrossRegionNonL4 = lower_triangle & ~pair_hasl4 & ~same_region;

% region-specific
regs_of_interest = ["V1", "M1", "V1L", "V1R"];
for kR = 1:length(regs_of_interest)
    reg = regs_of_interest(kR);
    in_region = contains(regions, reg) & contains(regions, reg)';
    
    s_masks.(sprintf('In%s', reg)) = ...
        lower_triangle & same_region & in_region;
    s_masks.(sprintf('In%sL4', reg)) = ...
        lower_triangle & pair_hasl4 & same_region & in_region;
    s_masks.(sprintf('In%sNonL4', reg)) = ...
        lower_triangle & ~pair_hasl4 & same_region & in_region;
end

end
