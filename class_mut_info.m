function minfo = class_mut_info(X)
% Similar to class_cond_entropy, but computes
% mututal information rather than conditional entropy.

cent = class_cond_entropy(X);

[n, m] = size(X);

% get full entropy of each class
ent = zeros(1, m);
all_classes = unique(X(:));

for k = all_classes(:)'
    pK = sum(X == k) / n;
    new_entropies = ent - pK .* log2(pK);
    
    % avoid nans
    bUpdate = pK > 0;
    ent(bUpdate) = new_entropies(bUpdate);
end

minfo = ent - cent;
    
end