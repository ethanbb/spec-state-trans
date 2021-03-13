function b_include_snippits = get_snippits_to_include(start_times, n_snippits, before, l, artifacts)
% Given the snippit start times and duration parameters and a matrix of artifact times,
% return a logical vector of which snippits should be kept (do not intersect with artifacts).

% first, trim the start times if necessary
real_starts = start_times - before;
real_starts(real_starts < 0) = [];
real_starts = real_starts(1:n_snippits);
ends = real_starts + l;

% algorithm: check how many artifacts each snippet is completely after and how many it is 
% completely before.
% snippet does not intersect an artifact iff the sum of these equals the number of artifacts.
n_artifacts = size(artifacts, 1);
n_after = sum(real_starts(:)' > artifacts(:, 2), 1);
n_before = sum(ends(:)' < artifacts(:, 1), 1);
b_include_snippits = (n_after + n_before == n_artifacts);

end
