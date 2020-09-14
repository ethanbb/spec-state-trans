function [pc_data_all, fh] = mt_pca(result_files, options)
% Do PCA in the frequency domain on multitaper results
% result_files: cell of mt_res mat files to use
% options: overrides of opts struct below
%
% the threshold types are:
%   'cumvar' - keep components until we reach this % cumulative explained
%              variance
%   'var'    - keep each component that explains at least this % variance
%   'comps'  - keep this number of (top) components
%
% returns a cell containing PC-space data for each input dataset
% if input is not a cell, output is unpacked one layer (but is still a cell
% of data from each channel)
% second output is an array of figure handles.

opts = struct(...
    'pxx_name',         'pxx',      ... variable name of input time-frequency data
    'name',             'pxx_pca',  ... variable name of output in each input file
    'chans',            'all',      ... of processed channels, which ones to use for each dataset (array of indices or 'all')
                                    ... can also be a cell of selections for each input
    'thresh_type',      'cumvar',   ... type of threshold for components to keep (see above)
    'thresh',           75,         ... threshold, interpreted according to thresh_type
    'save',             true        ... if true, saves 'pxx_pca' for each dataset within same file
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

if ~iscell(opts.chans)
    opts.chans = repmat({opts.chans}, size(result_files));
end

%------------ Load results and combine --------------%

unpack = false;
if ~iscell(result_files)
    unpack = true;
    result_files = {result_files};
end

n_files = numel(result_files);
c_pxx = cell(1, n_files);
c_lengths = cell(1, n_files);

nchans_orig = zeros(n_files, 1);

for kF = 1:n_files
    mfile = matfile(result_files{kF});
    nchans_orig(kF) = numel(mfile.(opts.pxx_name));
    
    if kF == 1
        % get the frequency axis while we're here
        % assume it's the same for all datasets, or else what are we
        % doing?
        
        freq_axis = mfile.freq_grid;
    end
    
    chans = opts.chans{kF};
    if ischar(chans)
        assert(strcmpi(chans, 'all'), 'Unrecognized ''chans'' input');
        opts.chans{kF} = 1:nchans_orig(kF);
        chans = opts.chans{kF};
    end
    
    pxx_file = mfile.(opts.pxx_name);
    c_pxx{kF} = [pxx_file{chans}];
    c_lengths{kF} = cellfun(@(p) size(p, 2), pxx_file(chans));
end
pxx = cell2mat(c_pxx);

%------------------ Do PCA ---------------------%

if strcmp(opts.thresh_type, 'comps')
    [coeff, pc_data, ~, ~, explained] = pca(pxx.', 'NumComponents', opts.thresh);
else
    [coeff, pc_data, ~, ~, explained] = pca(pxx.');
    
end

fh = gobjects(2, 1);

fh(1) = figure;
cumvar = cumsum(explained);
plot(cumvar, 'r*-', 'LineWidth', 2);
xlabel('PC#');
ylabel('Cumulative variance explained');
xlim([1, 50]);

if strcmp(opts.thresh_type, 'var')
    %-------- Apply variance cutoff ---------%
    
    num2keep = find(explained >= opts.thresh, 1, 'last');
    if num2keep == 0
        error('No components passed the explained variance threshold!');
    end
    
elseif strcmp(opts.thresh_type, 'cumvar')
    num2keep = find(cumvar > opts.thresh, 1);
else
    num2keep = opts.thresh;
end

%---------- Visualize the PCs ----------%

fh(2) = figure;
imagesc(coeff(:, 1:num2keep));
axis tight;
set(gca, 'YDir', 'normal');
set(gca, 'XTick', 1:num2keep);
xlabel('PC#');

if size(coeff, 1) == length(freq_axis) % otherwise don't bother with the y axis
    ytickinds = 1:20:size(coeff, 1);
    set(gca, 'YTick', ytickinds);
    set(gca, 'YTickLabel', freq_axis(ytickinds));
    ylabel('Frequency (Hz)');
end

if strcmp(opts.thresh_type, 'var')
    title(sprintf('Loadings of PCs each explaining >= %.2f%% of variance', opts.thresh));
elseif strcmp(opts.thresh_type, 'cumvar')
    title(sprintf('Loadings of PCs cumulatively explaining %.2f%% of variance', cumvar(num2keep)));
else
    title(sprintf('Loadings of top %d PCs', num2keep));
end

%------------- Reassemble cell of results -----------------%

n_assigned = 0;
pc_data_all = cell(size(result_files));

for kF = 1:n_files
    pc_data_all{kF} = cell(size(c_lengths{kF}));
    
    for kC = 1:numel(c_lengths{kF})
        n_pts = c_lengths{kF}(kC);
        pc_data_all{kF}{kC} = pc_data(n_assigned + (1:n_pts), 1:num2keep).';
        n_assigned = n_assigned + n_pts;
    end
end
assert(n_assigned == size(pc_data, 1), 'Uh oh, there''s an indexing problem');

%-------------- Save representations in individual files ------------%

if opts.save
    for kF = 1:n_files
        mfile = matfile(result_files{kF}, 'Writable', true');
        mfile.(opts.name) = cell(nchans_orig(kF), 1);
        mfile.(opts.name)(opts.chans{kF}, 1) = pc_data_all{kF};
        mfile.([opts.name, '_coeff']) = coeff(:, 1:num2keep);
        mfile.([opts.name, '_opts']) = opts;
    end
end

if unpack
    pc_data_all = pc_data_all{1};
end

end
