function [new_U, varargout] = rescale_to_pmf(U, varargin)
% Makes V or cell of Vs (from NMF) into a series of probability mass functions by the following
% steps:
%  * Divide each V by the maximum sum over components over all Vs (such that they each sum to <= 1
%    at each point), and multiply U by the same factor.
%  * Add a column to U consisting of all zeros
%  * Add a column to each V consisting of 1 - sum(V, 2) (such that all rows sum to 1).

scale_factor = max(cellfun(@(v) max(sum(v, 2)), varargin));
scale_factor = scale_factor + eps(scale_factor); % make sure there's no discrepancy when dividing each element
Vs = cellfun(@(v) v / scale_factor, varargin, 'uni', false);
U = U * scale_factor;

new_U = [U, zeros(size(U, 1), 1)];
varargout = cellfun(@(v) [v, 1 - sum(v, 2)], Vs, 'uni', false);

end
