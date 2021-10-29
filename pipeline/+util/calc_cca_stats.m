function [meanrho, redundancy_X, redundancy_Y] = calc_cca_stats(X, Y)
% Compute the mean rho (canoncial correlation) and global redundancy indices (asymmetric)
% between given N x P variable X and N x Q variable Y.

[~, ~, rhos, U, V] = canoncorr(X, Y);

meanrho = mean(rhos);

% Ref: https://github.com/andersonwinkler/toolbox/blob/master/share/redundancy.m
% Canonical loadings
At = corr(X, U);
Bt = corr(Y, V);

% Variance explained by each canonical variable
explU = mean(At.^2);
explV = mean(Bt.^2);

% Redundancy indices
rsq = rhos .^ 2;
Ru = explU .* rsq;
Rv = explV .* rsq;

% Global redundancy
redundancy_X = sum(Ru); % proportion of X variance explained
redundancy_Y = sum(Rv); % proportion of Y variance explained

end

