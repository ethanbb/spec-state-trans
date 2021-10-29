function [mean_err, err_each] = get_mean_nmf_recon_error(nmf_mfiles)
% Given a cell of NMF result matfile objects, get the mean
% NMF reconstruction error over all days and channels.

if ~isa(nmf_mfiles{1}, 'matlab.io.MatFile')
    nmf_mfiles = cellfun(@matfile, nmf_mfiles, 'uni', false);
end

n_chans_each = cellfun(@(mf) size(mf, 'all_chans', 1), nmf_mfiles(:));
chans_offset = [0; cumsum(n_chans_each(1:end-1))];
err_each = zeros(sum(n_chans_each), 1);

for kF = 1:length(nmf_mfiles)
    inds = chans_offset(kF) + (1:n_chans_each(kF));
    err_each(inds) = recon_err_from_mfile(nmf_mfiles{kF});
end
mean_err = mean(err_each);
end

function err_each = recon_err_from_mfile(nmf_mfile)
data = nmf_mfile.pxx_cat;
U = nmf_mfile.nmf_U;
U = U{1};
V = nmf_mfile.nmf_V;
V = V{1};

err_each = cellfun(@recon_err_one, data, U, V);
end

function err = recon_err_one(data, U, V)
err = norm(data - U*V', 'fro') / norm(data, 'fro');
% err = mean((data - U*V').^2, 'all') / var(data, 1, 'all');
end
