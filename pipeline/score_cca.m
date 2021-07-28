function score_cca(nmf_res_mfiles, channel_subsets, use_separate_runs, do_null)
% Use canonical correlation analysis to compare NMF scores from pairs of channels.
% Uses the mean of CCA correlation coefficients as the measure of similarity.
%
% nmf_res_mfiles is a cell of writable MatFile objects,
% including null model data, on which to operate.
%
% channel_subsets, if not empty, should be a struct. Each value is a set of channels to operate on -
% either a vector of indices into nmf_V/nmf_U or a cell of channel names, which can include 
% wildcards (* and ?) to select multiple channels. (A channel may appear in
% more than one such set). The corresponding key identifies the set and will be prepended to the
% saved variables in each matfile, along with the dist_type.
% The default is to use all channels and the name 'all'.
% 
% use_separate_runs is a boolean controlling whether to make a full matrix
% comparing each channel's scores in NMF run 1 to each channel in run 2 (if true), 
% or to make a half matrix (lower triangular) using just run 1 and comparing distinct channels.
% Default is false.
%
% do_null is a boolean for whether to make a "null model" version of the matrix. (Requires that
% gen_null_model_data has been run. If use_separate_runs is false, this means comparing shuffled
% score matrices to each other.)

arguments
    nmf_res_mfiles (:,1) cell
    channel_subsets (1,1) struct = struct('all', '*')
    use_separate_runs (1,1) logical = false
    do_null (1,1) logical = false
end

set_names = fieldnames(channel_subsets);
n_sets = length(set_names);
n_files = length(nmf_res_mfiles);
for kF = 1:n_files
    %% Get dataset info
    
    res_mfile = nmf_res_mfiles{kF};
    run_name = res_mfile.run_name;    
    chan_names = res_mfile.chan_names;
    Vs = res_mfile.nmf_V;
    Vis = Vs{1};
    if use_separate_runs
        Vjs = Vs{2};
    else
        Vjs = Vis;
    end
    
    if do_null
        Vs_null = res_mfile.null_V;
        Vi_null = Vs_null{1};
    end
    
    %% Loop over channel sets
    for kS = 1:n_sets
        % get actual channels depending on the format
        set_name = set_names{kS};
        set_val = channel_subsets.(set_name);
        if ischar(set_val) || isstring(set_val)
            set_val = {set_val};
        end
        
        if isnumeric(set_val)
            set_chans = set_val;
        elseif iscell(set_val)
            % take channels based on names
            regexps = regexptranslate('wildcard', set_val);
            joined_regexp = strjoin(regexps, '|');
            
            set_chans = find(~cellfun('isempty', regexp(chan_names, joined_regexp)));
        else
            error('Unrecognized channel subset format');
        end
        
        n_chans = length(set_chans);
        
        Vis_set = Vis(set_chans);
        Vjs_set = Vjs(set_chans);
        if do_null
            Vi_null_set = Vi_null(set_chans);
        end
        set_chan_names = chan_names(set_chans);
        
        cca_sim = nan(n_chans, n_chans);
        if do_null
            cca_sim_null = zeros(n_chans, n_chans);
        end
        
        for iC = 1:n_chans
            Vi_real = Vis_set{iC};
            if do_null
                Vi_null = Vi_null_set{iC};
            end
            
            if use_separate_runs
                n_j = n_chans;
            else
                n_j = iC - 1;
            end
            
            for jC = 1:n_j
                Vj = Vjs_set{jC};
                [~, ~, ps] = canoncorr(Vi_real, Vj);
                cca_sim(iC, jC) = mean(ps);
                
                if do_null
                    if use_separate_runs
                        [~, ~, ps] = canoncorr(Vi_null, Vj);
                    else
                        [~, ~, ps] = canoncorr(Vi_null, Vi_null_set{jC});
                    end
                    cca_sim_null(iC, jC) = mean(ps);
                end
            end
        end
        
        %% Save
        res_mfile.([set_name, '_cca_sim']) = cca_sim;
        if do_null
            res_mfile.([set_name, '_cca_sim_null']) = cca_sim_null;
        end
        res_mfile.([set_name, '_chans']) = set_chan_names;

    end
end

end
