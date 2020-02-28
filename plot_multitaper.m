function plot_multitaper(mt_res, bsave, chans)
% Input: output struct of multitaper_analysis (contents of mt_res.mat)
% If provided, 'bsave' specifies whether to save the figures in the current
%               directory.
% If provided, 'chans' specifies which among the available channels to use
%               (i.e. [1, 3] would be the first and third *analyzed*
%               channel.) Must be a vector of length 1 or 2.

n_chans_in = length(mt_res.options.chans);
assert(n_chans_in >= 1, 'No channels in input data');

if nargin < 3 || isempty(chans)
    chans = 1:n_chans_in;
else
    assert(max(chans) <= n_chans_in, 'Channels beyond number available requested');
    chans = chans(:).';
end
n_chans = length(chans);
    
if nargin < 2 || isempty(bsave)
    bsave = true;
end

h_fig = figure;
h_ax = gobjects(n_chans, 2);

chan_names = mt_res.options.chan_names;

for kC = 1:n_chans
    
    chan = chans(kC);
    
    % Plot power and normalized power spectra    
    pxx_db = 10*log10(mt_res.pxx{chan});
        
    h_ax(kC, 1) = subplot(n_chans, 2, kC*2 - 1);
    newplot;
    surface(mt_res.time_grid, mt_res.freq_grid, pxx_db, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log');
    axis tight;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Power in %s (dB)', chan_names{kC}));

    % include normalized/centered plot
    pxx_norm = mt_res.pxx{chan} ./ sum(mt_res.pxx{chan});
    pxx_norm_db = 10*log10(pxx_norm);
    pxx_norm_db_centered = pxx_norm_db - nanmean(pxx_norm_db, 2);
    
    h_ax(kC, 2) = subplot(n_chans, 2, kC*2);
    newplot;
    surface(mt_res.time_grid, mt_res.freq_grid, pxx_norm_db_centered, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log');
    axis tight;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Change from average of normalized spectrum in %s (dB)', chan_names{kC}));
end

linkaxes(h_ax, 'x');

% Save figure
if bsave
    savefig(h_fig, 'multitaper.fig', 'compact');
end

end