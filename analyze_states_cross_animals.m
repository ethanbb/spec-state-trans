% Response to JNeurosci reviewer 2 question about consistency of states across recording
% sessions/animals. Two parts:
%  * cluster data combined from all animals, for a set of channels, and see how much each animal
%    contributes to each cluster (use same NMF method as for main text)
%  * plot some characteristic data from each global state (depending on result of first part)

% Channels to compare across animals
% silly thing to deal with the fact that the same region has a different name in different recordings
channels = struct("M1V1", [
        struct('tag', 'M1_L4', 'chan_name', 'M1_L4')
        struct('tag', 'M1_Inf3', 'chan_name', 'M1_Inf3')
        struct('tag', 'V1_L4', 'chan_name', 'V1_L4')
        struct('tag', 'V1_Sup3', 'chan_name', 'V1_Sup3')
    ], "BilatV1", [
        struct('tag', 'V1_L4', 'chan_name', 'V1R_L4')
        struct('tag', 'V1_Sup3', 'chan_name', 'V1R_Sup3')
    ]);

sr_dirs = prepSR;
cd(sr_dirs.script);
inds_filename = 'characteristic_inds.mat';

if ~exist(inds_filename, 'file')
    disp([inds_filename ' not found; running script to generate']);
    get_characteristic_class_inds;
end

inds_mfile = matfile(inds_filename);
days_by_type = struct("M1V1", categorical(inds_mfile.m1v1_dates), 'BilatV1', categorical(inds_mfile.bilatv1_dates));
alldays = inds_mfile.run_dates;
topinds = inds_mfile.top_100_inds;

%% Load and combine data
% We end up with a big table

run_types = string(fieldnames(channels));
all_data = table(cell(length(run_types), 1), run_types, 'VariableNames', {'data', 'run_type'});

for kT = 1:length(run_types)
    this_chans = channels.(run_types(kT));
    n_chans = length(this_chans);
    chan_ids = categorical({this_chans.tag});
    
    this_days = days_by_type.(run_types(kT));
    n_days = length(this_days);
    all_data.data{kT} = table(cell(n_days, 1), this_days, 'VariableNames', {'data', 'day'});
    
    for kD = 1:n_days
        this_day = char(this_days(kD));
        day_ind = find(strcmp(this_day, alldays));
        topinds_day = topinds{day_ind};
        nmf_mfile = matfile(fullfile(this_day, 'nmf_res.mat'));
        chan_names = nmf_mfile.chan_names;
        chan_inds = cellfun(@(n) find(strcmp(n, chan_names)), {this_chans.chan_name});
        all_data.data{kT}.data{kD} = table(cell(n_chans, 1), chan_ids(:), 'VariableNames', {'data', 'channel'});
        
        % load the actual data
        for kC = 1:n_chans
            pxx_cat_1chan = nmf_mfile.pxx_cat(chan_inds(kC), 1);
            pxx_cat_1chan = pxx_cat_1chan{1};
            
            topinds_chan = topinds_day{kC};
            classes = 1:length(topinds_chan);
            classes(cellfun('isempty', topinds_chan)) = [];
            n_classes = length(classes);
            all_data.data{kT}.data{kD}.data{kC} = table(cell(n_classes, 1), classes', 'VariableNames', {'data', 'class'});
            
            for kK = 1:n_classes
                this_class = classes(kK);
                inds = topinds_chan{this_class};
                all_data.data{kT}.data{kD}.data{kC}.data{kK} = table(pxx_cat_1chan(:, inds)', 'VariableNames', {'data'});
            end
            
            all_data.data{kT}.data{kD}.data{kC} = explode_nested_tables(all_data.data{kT}.data{kD}.data{kC});
        end
        all_data.data{kT}.data{kD} = explode_nested_tables(all_data.data{kT}.data{kD});
    end
    all_data.data{kT} = explode_nested_tables(all_data.data{kT});
end
all_data = explode_nested_tables(all_data);


%% Now for each channel of interest, do NMF on whole concatenated dataset.
res_mfile = matfile('cross_animal_states.mat', 'Writable', true);
xval_fig_dir = fullfile('res_figs', 'cross_animal_nmf_xval');
cois = categories(all_data.channel);
res_mfile.channels = cois;
n_coi = length(cois);
n_comps = nan(n_coi, 1);

% break apart by channel
each_chan_data = cellfun(@(c) all_data(all_data.channel == c, :), cois, 'uni', false);

Us = cell(n_coi, 1);

parfor kC = 1:n_coi
    this_data = each_chan_data{kC};
    [this_n_comps, hfig] = util.nmf_ncomps_xval(this_data.data, 1, 15, 0.01);
    
    figure(hfig);
    title(sprintf('NMF cross-validation (%s)', cois{kC}));
    savefig(hfig, fullfile(xval_fig_dir, ['nmf_xval_', cois{kC}, '.fig']));
    n_comps(kC) = this_n_comps;
    close(hfig);
    
    % do NMF with inferred # of components
    [V, U] = sp_nnmf(this_data.data, this_n_comps, [], [], 500000);
    
    % do sorting and normalization as in concat_and_nmf (should probably be encapuslated)
    [~, peak_freqinds] = max(U);
    [~, order] = sort(peak_freqinds);
    U = U(:, order);
    V = V(:, order);
    
    % normalize
    norm_factor = vecnorm(U);
    U = U ./ norm_factor;
    V = V .* norm_factor;
    
    % get most likely "class"
    [~, classes] = max(V, [], 2);
    
    Us{kC} = U;
    this_data.nmf_V = V;
    this_data.global_class = classes;
    each_chan_data{kC} = this_data;
end

res_mfile.chan_names = cois;
res_mfile.n_comps = n_comps;
res_mfile.nmf_Us = Us;
res_mfile.data_by_chan = each_chan_data;
