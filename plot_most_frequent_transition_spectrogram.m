function mean_trans = plot_most_frequent_transition_spectrogram(rec_date, channel)
% Plot an average of the spectrograms 30 seconds before to 30 seconds after
% a state transition, picking the most frequent transition type.

% Load transition data to find the most frequent type
[~, input_s_all] = gather_exp_info;
this_input_s = input_s_all(strcmp(rec_date, {input_s_all.name}));
nmf_mfile = matfile(this_input_s.nmf_res_out);

all_chans = nmf_mfile.all_chans;
chan_ind = find(strcmp(all_chans, channel));
trans = nmf_mfile.filtered_transitions;
trans = trans{1}{chan_ind};

% zero out the diagonal
trans = trans .* (1 - eye(size(trans, 1)));

% find max row and column
[n_trans, ind] = max(trans, [], 'all', 'linear');
[from_state, to_state] = ind2sub(size(trans), ind);
fprintf('Using %d => %d (%d transitions)\n', from_state, to_state, n_trans);

% now make transition table to find the actual data
trans_table = get_state_transitions(this_input_s.mt_res_in, nmf_mfile, struct('save_filtered_classes', false));
btrans = trans_table.start_class == from_state & trans_table.end_class == to_state & trans_table.chan == chan_ind;
trans_table = trans_table(btrans, :);
assert(height(trans_table) == n_trans, 'Mismatch in # of transitions');

[~, ~, time_segs] = util.segmentize(this_input_s.mt_res_in, nmf_mfile.time_axis);

mt_fs = 10;
time_pre_post = 30;
trans_time_axis = -time_pre_post:1/mt_fs:time_pre_post;
n_time = length(trans_time_axis);
time_axis = nmf_mfile.time_axis;
freq_axis = nmf_mfile.freq_axis;

pxx = nmf_mfile.pxx_cat(chan_ind, 1);
pxx = pxx{1};
all_trans_spectra = nan(length(freq_axis), n_time, n_trans);

for kT = 1:n_trans
    kSeg = trans_table.segment(kT);
    trans_time = trans_table.time(kT);
    pre_time = trans_time - time_segs{kSeg}(1);
    post_time = time_segs{kSeg}(end) - trans_time;
    
    if pre_time < time_pre_post || post_time < time_pre_post
        fprintf('Skipping transition %d - not enough time before or after\n', kT);
        continue;
    end
    
    btime = time_axis >= trans_time - time_pre_post & time_axis <= trans_time + time_pre_post;
    all_trans_spectra(:, :, kT) = pxx(:, btime);
end

mean_trans = nanmean(all_trans_spectra, 3);

figure;
sanePColor(trans_time_axis, freq_axis, mean_trans, false, true);
set(gca, 'YScale', 'log');
xlabel('Time from transition (s)');
ylabel('Frequency (Hz)');
cb = colorbar;
cb.Label.String = 'Power rank';
title(sprintf('State %d to %d transitions', from_state, to_state));
hold on;
plot([0, 0], get(gca, 'YLim'), 'k--');

end
