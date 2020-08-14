function bicval = hmm_bic(model, data)
% Compute BIC of an HMM from hmmFit given the data it was fit on

logprob = hmmLogprob(model, data);
n = size(data, 2);

% get number of parameters depending on type of fit
k = numel(model.pi) + numel(model.A);
switch model.type
    case 'gauss'
        k = k + numel(model.emission.mu) + numel(model.emission.Sigma);
    case 'mixGaussTied'
        k = k + numel(model.emission.mu) + numel(model.emission.Sigma) + numel(model.emission.M);
    case 'discrete'
        k = k + numel(model.emission.T);
    case 'student'
        k = k + numel(model.emission.mu) + numel(model.emission.Sigma) + numel(model.emission.dof);
end

bicval = k * log(n) - 2 * logprob;

end