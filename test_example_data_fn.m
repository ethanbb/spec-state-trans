function test_example_data_fn(rec_date, channel, mt_file_paths, lfp_inds)
% Check whether the segments found by
% plot_example_data_from_classes are actually from the intended classes.

Fs = 1000;

[~, input_s_all] = gather_exp_info;
this_input_s = input_s_all(strcmp(rec_date, {input_s_all.name}));
nmf_mfile = matfile(this_input_s.nmf_res_out);
bChan = strcmp(channel, nmf_mfile.chan_names);
classes_all = nmf_mfile.nmf_classes;
classes_all = classes_all{1};
chan_classes = classes_all{bChan};

% construct a vector of selected class inds
n_files = length(this_input_s.mt_res_in);
taken_class = cell(n_files, 1);
seg_startsamples = cell(n_files, 1);
seg_endsamples = cell(n_files, 1);

for kF = 1:n_files
    mt_mfile = matfile(this_input_s.mt_res_in{kF});
    opts = mt_mfile.options;
    step_samps = opts.winstep * Fs;
    
    seg_samples = mt_mfile.seg_samples;
    seg_startsamples{kF} = cellfun(@(segsamps) segsamps(1), seg_samples);
    seg_endsamples{kF} = cellfun(@(segsamps) segsamps(end), seg_samples);
    seg_nwindows = cellfun('length', mt_mfile.seg_windows);
%     seg_cum_nwindows{kF} = cumsum(seg_nwindows);
    taken_class{kF} = arrayfun(@(nsamp) zeros(nsamp, 1), seg_nwindows(:), 'uni', false);
end

[n_samples, n_classes] = size(mt_file_paths);
for kK = 1:n_classes
    for kS = 1:n_samples
        bfile = strcmp(mt_file_paths{kS, kK}, this_input_s.mt_res_in);
        start_samp = lfp_inds(1, kS, kK);
        end_samp = lfp_inds(2, kS, kK);
        
        % locate the windows
        kseg = find(start_samp >= seg_startsamples{bfile}, 1, 'last');
        assert(find(end_samp <= seg_endsamples{bfile}, 1) == kseg, ...
            'Hmm, sample crosses segment boundary');
        
        kwin_first = floor((start_samp - seg_startsamples{bfile}(kseg)) / step_samps);
        kwin_last = floor((seg_endsamples{bfile}(kseg) - end_samp) / step_samps);
        taken_class{bfile}{kseg}(1+kwin_first:end-kwin_last) = kK;
    end
end

% concatenate, then compare to saved classes
taken_classes_perseg = vertcat(taken_class{:});
taken_classes_all = vertcat(taken_classes_perseg{:});

assert(length(taken_classes_all) == length(chan_classes), 'Hmm, class vectors are of different lengths');
btaken = taken_classes_all ~= 0;

if all(taken_classes_all(btaken) == chan_classes(btaken))
    disp('All good, the classes match');
else
    error('Class mismatch - investigate.');
end

end
