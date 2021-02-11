function plot_multitaper(result, options)
% Input: output struct of multitaper_analysis (contents of mt_res.mat),
% MatFile object, or MAT-filename containing these results.
% If a cell is passed, concatenates them together across time.

default_savedir = pwd;
if ~iscell(result)
    result = {result};
end
result = result(:);

if ischar(result{1})
    default_savedir = fileparts(result{1});
    result = cellfun(@matfile, result, 'uni', false);
end

opts = struct(...
    'pxx_name',  'pxx_pp',               ... fieldname (within result) of data to plot, or cell to plot multiple
    'take_log',  false,                  ... logical array same size as pxx_name of whether to take log first
    'clim',      [],                     ... color limits (default = automatic)
    'label',     [],                     ... label for each field (entry of pxx_name)
    'chans',     'all',                  ... which channels, of those analyzed, to plot
    'save',      false,                  ... whether to save the figure
    'savedir',   default_savedir,        ... directory, if saving
    'filename',  'multitaper.fig',       ... filename, if saving
    'remove_nans', false,                ... if true, take out nan segments and plot on an x-axis of just valid data
    'axes',      []                      ... if provided, plots into these axes instead of making a new figure
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

% input validation
n_inputs = length(result);
for kR = 1:n_inputs
    opts = validate_res_and_opts(result{kR}, opts);
end
n_chans = length(opts.chans);
n_fields = length(opts.pxx_name);
if ~isempty(opts.axes)
    assert(numel(opts.axes) == n_chans * n_fields, 'Wrong # of axes provided');
end

% get data
pxx = repmat({cell(1, n_inputs)}, n_chans, n_fields);
time_s = 0;
for kI = 1:n_inputs
    this_res = result{kI};
    for kF = 1:n_fields
        this_field = this_res.(opts.pxx_name{kF});
        time_grid = this_res.time_grid;
        if ~exist('dt', 'var')
            dt = mode(diff(time_grid));
        else
            assert(mode(diff(time_grid)) - dt < 1e-10, 'Inputs have different sample rates');
        end

        for kC = 1:n_chans
            chan_data = this_field{opts.chans(kC)};
            if kF == 1 && kC == 1
                if opts.remove_nans
                    this_chan_keep = ~all(isnan(chan_data));
                else
                    this_chan_keep = true(1, size(chan_data, 2));
                end
                time_s = time_s + sum(this_chan_keep) * dt;
            end
            
            pxx{kC, kF}{kI} = chan_data(:, this_chan_keep);
        end
    end
end
pxx = cellfun(@cell2mat, pxx, 'uni', false);


% set up figure
if isempty(opts.axes)
    h_fig = figure('Position', [0, 0, 900*n_fields, 450*n_chans]);
    h_ax = gobjects(n_chans, n_fields);
    for kC = 1:n_chans
        for kF = 1:n_fields
            h_ax(kC, kF) = subplot(n_chans, n_fields, (kC-1)*n_fields + kF);
        end
    end
else
    h_ax = reshape(opts.axes, n_chans, n_fields);
end

for kC = 1:n_chans
    for kF = 1:n_fields
        axes(h_ax(kC, kF)); % have to do this b/c sanePColor doesn't support Axes input

        this_pxx = pxx{kC, kF};
        if opts.take_log(kF)
            this_pxx = 10*log10(this_pxx);
        end

        sanePColor(dt:dt:time_s, opts.freq_grid, this_pxx, false, true);
        set(gca, 'YScale', 'log');
        set(gca, 'Interactions', regionZoomInteraction); % dataTipInteraction causes lots of lag
        if ~isempty(opts.clim)
            set(gca, 'CLim', opts.clim);
        end
        axis tight;
        if opts.remove_nans
            xlabel('Time (excluding artifacts) (s)');
        else
            xlabel('Time (s)');
        end
        ylabel('Frequency (Hz)');
        title_str = sprintf('%s in %s', opts.label{kF}, opts.chan_names{kC});
        if opts.take_log(kF)
            title_str = sprintf('Log of %s', title_str);
        end
        title(title_str, 'Interpreter', 'none');
    end
end

linkaxes(h_ax, 'x');

% Save figure
if opts.save
    if ~exist('h_fig', 'var')
        warning('Did not save figure since axes were provided');
        return;
    end
    if ~exist(opts.savedir, 'dir')
        % try to create it
        if ~mkdir(opts.savedir)
            error('Could not create save directory %s', opts.savedir);
        end
    end

    savefig(h_fig, fullfile(opts.savedir, opts.filename), 'compact');
end

end

function opts = validate_res_and_opts(result, opts)
% Here, result is either a struct or a MatFile object.

all_chan_names = result.name;
n_chans_in = length(all_chan_names);
assert(n_chans_in >= 1, 'No channels in input data');

% process options
if ischar(opts.chans)
    assert(strcmpi(opts.chans, 'all'), 'Unrecognized ''chans'' input');
    opts.b_chans_all = true;
    opts.chans = 1:n_chans_in;
elseif isfield(opts, 'b_chans_all') && opts.b_chans_all
    assert(all(opts.chans == 1:n_chans_in), 'Inputs have different #s of channels');
else
    assert(max(opts.chans) <= n_chans_in, 'Channels beyond number available requested');
end
opts.chans = opts.chans(:).';

chan_names = all_chan_names(opts.chans);
if ~isfield(opts, 'chan_names')
    opts.chan_names = chan_names;
else
    assert(length(opts.chan_names) == length(chan_names) && all(strcmp(opts.chan_names, chan_names)), ...
        'Channel names don''t match between inputs');
end

freq_grid = result.freq_grid;
if ~isfield(opts, 'freq_grid')
    opts.freq_grid = freq_grid;
else
    assert(length(opts.freq_grid) == length(freq_grid) && all(opts.freq_grid == freq_grid), ...
        'Frequencies don''t match between inputs');
end

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
% use some smart default labels
predef_labels = struct(...
    'pxx',          'Power', ...
    'phase',        'Phase', ...
    'pxx_pp',       'Preprocessed power', ...
    'pxx_smooth',   'Smoothed power', ...
    'pxx_rankord',  'Rank-order power');

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
end
