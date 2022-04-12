function test_mt_timeres
% 2-point discrimination test for time resolution of multitaper analysis.
% Use util.sim_mt_and_smoothing on sample data with 2 impulses to determine how close they have to
% be for their peaks to be indistinguishable.

start_gap_s = 6; % seconds; this should be a gap where they're guaranteed to be distinguishable
Fs = 1000;
possible_gaps = 1:start_gap_s*Fs; % samples
gaps_tested = false(1, length(possible_gaps));
sig_length = 180; % seconds
sig_time = -sig_length/2:1/Fs:sig_length/2;
midpt_samp = find(sig_time == 0);

gap = possible_gaps(end);
while true
    fprintf('Trying gap: %g s\n', gap / Fs);
    data = make_sample_data(gap);
    spectrogram = util.sim_mt_and_smoothing(data);
    b_dist = is_distinguishable(spectrogram);
    gaps_tested(gap) = true;
    
    if b_dist  % distinguishable; decrease the gap
        lower_gap = find(gaps_tested(1:gap-1), 1, 'last');
        if isempty(lower_gap)
            lower_gap = 0;
        end
        
        gap_diff = gap - lower_gap;
        if gap_diff == 1
            % we are at the minimum distinguishable threshold
            break;
        end
        
        % halve the distance between here and the next lower indistinguishable gap
        gap = gap - floor(gap_diff / 2);
        
    else  % indistinguishable
        upper_gap = gap + find(gaps_tested(gap+1:end), 1);
        if isempty(upper_gap)
            error('Start gap not high enough to be distinguishable');
        end
        
        gap_diff = upper_gap - gap;
        if gap_diff == 1
            % next up is the minimum distinguishable gap
            gap = gap + 1;
            break;
        end
        
        % halve the distance between there and the next higher distinguishable gap
        gap = gap + floor(gap_diff / 2);
    end
end

min_dist_gap_secs = gap / Fs;
fprintf('The minimum distinguishable gap is %g s.\n', min_dist_gap_secs);


    function data = make_sample_data(gap)
        data = zeros(1, length(sig_time));
        data(midpt_samp - floor(gap/2)) = 1;
        data(midpt_samp + ceil(gap/2)) = 1;
    end

    function b_distinguishable = is_distinguishable(spectrogram)
        n_max = sum(islocalmax(spectrogram, 2), 2);
        b_distinguishable = any(n_max > 1);
    end

end
