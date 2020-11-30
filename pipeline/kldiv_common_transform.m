function [X1, X2] = kldiv_common_transform(V1, V2)
% Input: cells of scores from each channel, from 2 different NMF runs (1 and 2)
% Output: cells of transform matrices (of same size as input) that minimize
%  D_KL(V1{i}X1{i}, V2{j}X2{j}) for each pair (i,j).

C1 = length(V1);
C2 = length(V2);

% Get number of components in each V
ks1 = cellfun('size', V1, 2);
ks2 = cellfun('size', V2, 2);

% for indexing into the full X
k1_offsets = cumsum([0; ks1(1:end)]);
k2_offsets = sum(ks1) + cumsum([0; ks2(1:end)]);

% Dimension of each score vector after transformation
k = max([ks1(:); ks2(:)]);

total_ks = sum(ks1) + sum(ks2);
% "X" for fmincon will be just each X1 and X2 stacked vertically.
% (size: total_ks x k)
% Equality constraints: each row of X should sum to 1, such that any affine
% combination of them also sums to 1.
Aeq = repmat(eye(total_ks), 1, k);
beq = ones(total_ks, 1);
x0 = rand(total_ks, k);
x0 = x0 ./ sum(x0, 2);

% Bounds on X - this may not be strictly necessary, but it's a lot more
% feasible than directly constraining each element of each VX and it
% guarantees the outputs are probabilities
lb = zeros(size(x0));
ub = ones(size(x0));

optimopts = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e7);

X = fmincon(@kldiv_common_objective, x0, [], [], Aeq, beq, lb, ub, [], optimopts);

% Now split apart outputs
Xs = mat2cell(X, [ks1; ks2]);
X1 = Xs(1:C1);
X2 = Xs(C1+1:end);

    function kldiv = kldiv_common_objective(X)
        kldiv = 0;
        for ci = 1:C1
            Vi = V1{ci};
            Xi = X(k1_offsets(ci) + (1:ks1(ci)), :);
            VXi = Vi * Xi;
            
            for cj = 1:C2
                Vj = V2{cj};
                Xj = X(k2_offsets(cj) + (1:ks2(cj)), :);
                VXj = Vj * Xj;
                
                kldiv = kldiv + sum(sum(VXi .* (log(VXi) - log(VXj))));
            end
        end
    end

end