function score_cca(nmf_res_mfiles, channel_subsets)
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

arguments
    nmf_res_mfiles (:,1) cell
    channel_subsets (1,1) struct = struct('all', '*')
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
    Vjs = Vs{2};
    Vs_null = res_mfile.null_V;
    Vi_null = Vs_null{1};
    
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
        Vi_null_set = Vi_null(set_chans);
        set_chan_names = chan_names(set_chans);
        
        cca_sim = zeros(n_chans, n_chans);
        cca_sim_null = zeros(n_chans, n_chans);
        
        for iC = 1:n_chans
            Vi_real = Vis_set{iC};
            Vi_null = Vi_null_set{iC};
            
            for jC = 1:n_chans
                Vj = Vjs_set{jC};
                [~, ~, ps] = canoncorr(Vi_real, Vj);
                cca_sim(iC, jC) = mean(ps);
                
                [~, ~, ps] = canoncorr(Vi_null, Vj);
                cca_sim_null(iC, jC) = mean(ps);
            end
        end
        
        %% Save
        res_mfile.([set_name, '_cca_sim']) = cca_sim;
        res_mfile.([set_name, '_cca_sim_null']) = cca_sim_null;
        res_mfile.([set_name, '_chans']) = set_chan_names;
        
        %% Plot matrix of similarity
        mfile_dir = fileparts(res_mfile.Properties.Source);
        
        hf = plot_dist_mat(cca_sim, set_chan_names, ...
            sprintf('%s - %s', run_name, set_name), 'cca_mean_r');
        savefig(hf, fullfile(mfile_dir, sprintf('cca_%s_%s.fig', run_name, set_name)));
        
        hf2 = plot_dist_mat(cca_sim_null, set_chan_names, ...
            sprintf('%s - %s (null model)', run_name, set_name), 'cca_mean_r');
        savefig(hf2, fullfile(mfile_dir, sprintf('cca_null_%s_%s.fig', run_name, set_name)));
    end
end

end
