function plot_multitaper(result, options)
% Input: output struct of multitaper_analysis (contents of mt_res.mat),
% MatFile object, or MAT-filename containing these results.

default_savedir = pwd;
if ischar(result)
    default_savedir = fileparts(result);
    result = load(result);
end

% use some smart default labels
predef_labels = struct(...
    'pxx',          'Power', ...
    'phase',        'Phase', ...
    'pxx_pp',       'Preprocessed power', ...
    'pxx_smooth',   'Smoothed power', ...
    'pxx_rankord',  'Rank-order power');

opts = struct(...
    'pxx_name',  'pxx_pp',               ... fieldname (within result) of data to plot, or cell to plot multiple
    'take_log',  false,                  ... logical array same size as pxx_name of whether to take log first
    'clim',      [],                     ... color limits (default = automatic)
    'label',     [],                     ... label for each field (entry of pxx_name)
    'chans',     'all',                  ... which channels, of those analyzed, to plot
    'save',      false,                  ... whether to save the figure
    'savedir',   default_savedir,        ... directory, if saving
    'filename',  'multitaper.fig'        ... filename, if saving
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

mt_opts = result.options;
n_chans_in = length(mt_opts.chans);
assert(n_chans_in >= 1, 'No channels in input data');

% process options
if ischar(opts.chans)
    assert(strcmpi(opts.chans, 'all'), 'Unrecognized ''chans'' input');
    opts.chans = 1:n_chans_in;
else
    assert(max(opts.chans) <= n_chans_in, 'Channels beyond number available requested');
end
opts.chans = opts.chans(:).';
chan_names = mt_opts.chan_names(opts.chans);
n_chans = length(opts.chans);

% validate pxx_name
if ischar(opts.pxx_name)
    opts.pxx_name = {opts.pxx_name};
end
opts.pxx_name = opts.pxx_name(:);
n_fields = length(opts.pxx_name);

if isstruct(result)
    has_all = all(isfield(result, opts.pxx_name));
else
    has_all = all(cellfun(@(f) isprop(result, f), opts.pxx_name));
end
assert(has_all, 'Field name(s) to plot not present in data');

% validate labels and 'take_log'
if ischar(opts.label)
    assert(n_fields == 1, 'Not enough field labels (expected %d)', n_fields);
    opts.label = {opts.label};

elseif isempty(opts.label)
    % try to assign reasonable defaults
    opts.label = cell(size(opts.pxx_name));
    for kF = 1:n_fields
        if isfield(predef_labels, opts.pxx_name{kF})
            opts.label{kF} = predef_labels.(opts.pxx_name{kF});
        else
            % default to just the field name
            opts.label{kF} = opts.pxx_name{kF};
        end
    end

else
    opts.label = opts.label(:);
    assert(length(opts.label) == n_fields, 'Wrong # of field labels (expected %d)', n_fields);
end

opts.take_log = opts.take_log(:);
if isscalar(opts.take_log)
    opts.take_log = repmat(opts.take_log, n_fields, 1);
end
assert(length(opts.take_log) == n_fields, 'Wrong # of ''take_log'' entries (expected %d)', n_fields);

% get data
pxx = cell(n_chans, n_fields);
for kF = 1:n_fields
    this_field = result.(opts.pxx_name{kF});

    for kC = 1:n_chans
        pxx{kC, kF} = this_field{opts.chans(kC)};
    end
end

% set up figure
h_fig = figure('Position', [0, 0, 900*n_fields, 450*n_chans]);
h_ax = gobjects(n_chans, n_fields);

for kC = 1:n_chans
    for kF = 1:n_fields
        h_ax(kC, kF) = subplot(n_chans, n_fields, (kC-1)*n_fields + kF);

        this_pxx = pxx{kC, kF};
        if opts.take_log(kF)
            this_pxx = 10*log10(this_pxx);
        end

        sanePColor(result.time_grid, result.freq_grid, this_pxx, false, true);
        set(gca, 'YScale', 'log');
        set(gca, 'Interactions', regionZoomInteraction); % dataTipInteraction causes lots of lag
        if ~isempty(opts.clim)
            set(gca, 'CLim', opts.clim);
        end
        axis tight;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(sprintf('%s in %s', opts.label{kF}, chan_names{kC}), 'Interpreter', 'none');
    end
end

linkaxes(h_ax, 'x');

% Save figure
if opts.save
    if ~exist(opts.savedir, 'dir')
        % try to create it
        if ~mkdir(opts.savedir)
            error('Could not create save directory %s', opts.savedir);
        end
    end

    savefig(h_fig, fullfile(opts.savedir, opts.filename), 'compact');
end

end
