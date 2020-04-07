function plot_multitaper(result, options)
% Input: output struct of multitaper_analysis (contents of mt_res.mat), or
% MAT-filename containing these results.

default_savedir = pwd;
if ischar(result)
    default_savedir = fileparts(result);
    result = load(result);
end

opts = struct(...
    'pxx_name',  'pxx',                  ... fieldname (within result) of data to plot
    'chans',     'all',                  ... which channels, of those analyzed, to plot
    'save',      true,                   ... whether to save the figure
    'savedir',   default_savedir,        ... directory, if saving
    'filename',  'multitaper.fig',       ... filename, if saving
    'normonly',  true,                   ... only make the normalized plot
    'xfreqnorm', false                   ... true to normalize pwr across freqs at each time
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

n_chans_in = length(result.options.chans);
assert(n_chans_in >= 1, 'No channels in input data');

% process options
if ischar(opts.chans)
    assert(strcmpi(opts.chans, 'all'), 'Unrecognized ''chans'' input');
    opts.chans = 1:n_chans_in;
else
    assert(max(opts.chans) <= n_chans_in, 'Channels beyond number available requested');
end
opts.chans = opts.chans(:).';
n_chans = length(opts.chans);

% get subset of pxx and names according to chans
pxx = result.(opts.pxx_name)(opts.chans);
chan_names = result.options.chan_names(opts.chans);

% set up figure
n_cols = 2 - opts.normonly;
h_fig = figure('Position', [0, 0, 900*n_cols, 450*n_chans]);
h_ax = gobjects(n_chans, n_cols);

for kC = 1:n_chans
    
    pxx_db = 10*log10(pxx{kC});

    if ~opts.normonly
        % Plot power and normalized power spectra

        h_ax(kC, 2) = subplot(n_chans, n_cols, kC*2 - 1);
        newplot;
        surface(result.time_grid, result.freq_grid, pxx_db, 'EdgeColor', 'none');
        set(gca, 'YScale', 'log');
        axis tight;
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title(sprintf('Power in %s (dB)', chan_names{kC}));
    end

    % include normalized/centered plot
    if opts.xfreqnorm
        pxx_norm = pxx{kC} ./ sum(pxx{kC});
        pxx_norm_db = 10*log10(pxx_norm);
        pxx_centered = pxx_norm_db - nanmean(pxx_norm_db, 2);
    else % just subtract mean in log space (i.e. normalize by geometric mean)
        pxx_centered = pxx_db - nanmean(pxx_db, 2);
    end
    
    h_ax(kC, 1) = subplot(n_chans, n_cols, kC*n_cols);
    newplot;
    surface(result.time_grid, result.freq_grid, pxx_centered, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log');
    axis tight;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    if opts.xfreqnorm
        fs = 'Normalized power distribution over time in %s (dB)';
    else
        fs = 'Normalized power in %s (dB)';
    end
    title(sprintf(fs, chan_names{kC}));
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