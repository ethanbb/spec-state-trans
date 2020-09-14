function gen_null_model_data(nmf_res_mfiles)
% Given NMF results, for each channel on each day, fit a Markov chain to the sequence of
% largest-score classes, generate a new sequence from this, and randomly sample score vectors with
% replacement from original data to make a null-model NxK score matrix.
%
% nmf_res_mfiles is a cell of writable MatFile objects on which to operate.

n_files = length(nmf_res_mfiles);
for kF = 1:n_files
    rstate = rng('shuffle');
    
    % Get data for all reps/channels
    mfile = nmf_res_mfiles{kF};
    mfile.null_rstate = rstate;
    real_V = mfile.nmf_V;
    null_V = real_V; % just to copy shape
    real_classes = mfile.nmf_classes;
    null_classes = real_classes; % ditto
    trans_allreps = mfile.nmf_transitions;
    
    n_times = length(mfile.time_axis);
    n_reps = length(trans_allreps);

    for kR = 1:n_reps
        trans_allchans = trans_allreps{kR};
        n_chans = length(trans_allchans);
        
        for kC = 1:n_chans
            trans = trans_allchans{kC};
            n_classes = size(trans, 1);
            
            % Filter out classes that are never present (i.e. components that are never the max)
            valid_classes = find(sum(trans, 2) > 0);
            trans_valid = trans(valid_classes, valid_classes);
            
            % Initial state
            start_state = 0;
            while ~ismember(start_state, valid_classes)
                start_state = randsample(real_classes{kR}{kC}, 1);
            end
            x0 = zeros(1, length(valid_classes));
            x0(valid_classes == start_state) = 1;
            
            % Make Markov chain based on transition data
            model = dtmc(trans_valid);
            null_class_inds = simulate(model, n_times-1, 'X0', x0);
            null_classes{kR}{kC} = valid_classes(null_class_inds);
            
            % Sample from real V to make null V
            for kK = 1:n_classes
                inds_byclass = find(real_classes{kR}{kC} == kK);
                null_inclass = null_classes{kR}{kC} == kK;
                null_V_inds = randsample(inds_byclass, sum(null_inclass), true);
                null_V{kR}{kC}(null_inclass, :) = real_V{kR}{kC}(null_V_inds, :);
            end          
        end
    end
    
    mfile.null_V = null_V;
    mfile.null_classes = null_classes;
end

end
