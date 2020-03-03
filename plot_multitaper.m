function plot_multitaper(result, options)
% Input: output struct of multitaper_analysis (contents of mt_res.mat), or
% MAT-filename containing these results.

default_savedir = pwd;
if ischar(result)
    default_savedir = fileparts(result);
    result = load(result);
end

opts = struct(...
    'chans',     'all',                  ... which channels, of those analyzed, to plot
    'save',      true,                   ... whether to save the figure
    'savedir',   default_savedir,        ... directory, if saving
    'filename',  'multitaper.fig'        ... filename, if saving
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{KO}) = options.(opts_in{kO});
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

h_fig = figure;
h_ax = gobjects(n_chans, 2);

chan_names = result.options.chan_names;

for kC = 1:n_chans
    
    chan = opts.chans(kC);
    
    % Plot power and normalized power spectra    
    pxx_db = 10*log10(result.pxx{chan});
        
    h_ax(kC, 1) = subplot(n_chans, 2, kC*2 - 1);
    newplot;
    surface(result.time_grid, result.freq_grid, pxx_db, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log');
    axis tight;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Power in %s (dB)', chan_names{kC}));

    % include normalized/centered plot
    pxx_norm = result.pxx{chan} ./ sum(result.pxx{chan});
    pxx_norm_db = 10*log10(pxx_norm);
    pxx_norm_db_centered = pxx_norm_db - nanmean(pxx_norm_db, 2);
    
    h_ax(kC, 2) = subplot(n_chans, 2, kC*2);
    newplot;
    surface(result.time_grid, result.freq_grid, pxx_norm_db_centered, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log');
    axis tight;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Change from average of normalized spectrum in %s (dB)', chan_names{kC}));
end

linkaxes(h_ax, 'x');

% Save figure
if opts.save
    savefig(h_fig, fullfile(opts.savedir, opts.filename), 'compact');
end

end