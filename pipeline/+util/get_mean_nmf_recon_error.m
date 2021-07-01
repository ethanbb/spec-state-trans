function mean_err = get_mean_nmf_recon_error(nmf_mfiles)
% Given a cell of NMF result matfile objects, get the mean
% NMF reconstruction error over all days and channels.

n_chans_total = sum(cellfun(@(mf) size(mf, 'all_chans', 1), nmf_mfiles));
err_sum = sum(cellfun(@recon_err_from_mfile, nmf_mfiles));
mean_err = err_sum / n_chans_total;

end

function err_sum = recon_err_from_mfile(nmf_mfile)
data = nmf_mfile.pxx_cat;
U = nmf_mfile.nmf_U;
U = U{1};
V = nmf_mfile.nmf_V;
V = V{1};

err_sum = sum(cellfun(@recon_err_one, data, U, V));
end

function err = recon_err_one(data, U, V)
err = norm(data - U*V', 'fro') / norm(data, 'fro');
end
