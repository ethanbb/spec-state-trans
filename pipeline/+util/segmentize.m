function [seg_starts, seg_ends, varargout] = segmentize(mt_res_files, varargin)
% Given a cell of multitaper results files (in order), divide each non-cell input into a cell
% of parts of lengths matching the (concatenated) clean semgents of multitaper windows.
% Each cell input will have this done to each element. Non-vector inputs will be divided
% along the first dimension.

seg_lengths = cell(1, length(mt_res_files));
for kR = 1:length(mt_res_files)
    mt_res = mt_res_files{kR};
    if ischar(mt_res)
        mt_res = matfile(mt_res);
    end
    seg_lengths{kR} = cellfun('length', mt_res.seg_windows);
end
seg_lengths = cell2mat(seg_lengths);
seg_ends = cumsum(seg_lengths);
seg_starts = [1, seg_ends(1:end-1) + 1];

varargout = segmentize_helper(varargin, seg_lengths);

end

function var_divided = segmentize_helper(var_in, seg_lengths)
if iscell(var_in)
    var_divided = cellfun(@(v) segmentize_helper(v, seg_lengths), var_in, 'uni', false);
else
    if isvector(var_in)
        var_in = var_in(:);
    end
    var_divided = mat2cell(var_in, seg_lengths);
end
end
