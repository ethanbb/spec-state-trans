function chan_fixup(mt_res, savedir)

merge_mt_res = @(old, new) cat(find(size(old) > 1, 1), old(1), new(1:2), old(3:4), new(3:4), old(6));

mt_res_old = matfile(fullfile(savedir, 'mt_res.mat'));
mt_res_old_opts = mt_res_old.options;

mt_res.options.chans = merge_mt_res(mt_res_old_opts.chans, mt_res.options.chans);
mt_res.options.chan_names = merge_mt_res(mt_res_old_opts.chan_names, mt_res.options.chan_names);
mt_res.name = merge_mt_res(mt_res_old.name, mt_res.name);
mt_res.pxx = merge_mt_res(mt_res_old.pxx, mt_res.pxx);

save(fullfile(savedir, 'mt_res.mat'), '-struct', 'mt_res', '-v7.3');

end