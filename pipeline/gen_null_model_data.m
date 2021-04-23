function gen_null_model_data(nmf_res_mfiles, reshuffle_rng)
% Given NMF results, for each channel on each day, fit a Markov chain to the sequence of
% largest-score classes, generate a new sequence from this, and randomly sample score vectors with
% replacement from original data to make a null-model NxK score matrix.
%
% nmf_res_mfiles is a cell of writable MatFile objects on which to operate.

n_files = length(nmf_res_mfiles);
for kF = 1:n_files
    mfile = nmf_res_mfiles{kF};

    if ~exist('reshuffle_rng', 'var') || isempty(reshuffle_rng) || reshuffle_rng ...
            || ~isprop(mfile, 'null_rstate')
        rng('shuffle');
        mfile.null_rstate = rng;
    else  % reuse the same random state
        rng(mfile.null_rstate);
    end
    
    % Get data for all reps/channels
    real_V = mfile.nmf_V;
    null_V = real_V; % just to copy shape
    real_classes = mfile.nmf_classes;
    null_classes = real_classes; % ditto
    transitions = mfile.nmf_transitions;
    
    for kR = 1:length(real_classes)
        for kC = 1:length(real_classes{kR})
            [null_V{kR}{kC}, null_classes{kR}{kC}] = util.shuffle_scores_markov(...
                real_V{kR}{kC}, real_classes{kR}{kC}, transitions{kR}{kC}, [], false);         
        end
    end
    
    mfile.null_V = null_V;
    mfile.null_classes = null_classes;
end

end
