function [pca_change, pca_change_time] = pca_change(result_file, options)
% For each mt_res.mat result file in the result_files input cell,
% for each selected channel (all processed ones by default),
% compute the magnitude of spectral change (i.e. change speed) in PCA
% space. By default, plots this along with a smoothed version, and also
% plots the cross-correlation between change speeds of the first included
% channels if there are at least 2. These defaults can be overriden by
% passing an options override struct as the second argument.
%
% First output is a nchans x time array of the change signal.
% Second output is the corresponding x axis (in seconds).

opts = struct(...
    'pca_name',      'pxx_pca',      ... variable name of PCA input within input file
    'rel_chans',     'all',          ... which channels to use, of those available (array or 'all' (default))
    'comps',         'all',          ... which components to use (array or 'all')
    'smooth_span',   10,             ... span in seconds to smooth the change (before interpolating)
    'smooth_method', 'movmean',      ... smoothing method (3rd argument of smoothdata)
    'norm_type',     2,              ... second input to vecnorm (default is Euclidian distance)
    'interp_factor', 1,              ... set to >1 to interpolate change by an integer factor
    'interp_method', 'makima',       ... interpolation method (4th argument of interp1)
    'xcorr_pairs',   {{}},           ... pairs (2-element vectors) of indices of rel_chans to plot xcorr of
    'plot',          true,           ... whether to do plotting
    'plot_raw',      false,          ... whether to include raw (unsmoothed) change in plot
    'save',          true,           ... whether to save result in original file
    'name',          'pca_change',   ... variable name to save (with '_time' appended for x-axis)
    'savefigs',      true            ... whether to save figures in parent dir of input file
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

%------------ validate options/input ------------------

data_s = load(result_file);

assert(isfield(data_s, opts.pca_name), 'Specified PCA results do not exist');

valid_chans = ~cellfun(@isempty, data_s.(opts.pca_name));
valid_names = data_s.name(valid_chans);
valid_data = data_s.(opts.pca_name)(valid_chans);

n_chans = numel(valid_data);
assert(n_chans > 0, 'Specified PCA results do not exist');

if ischar(opts.rel_chans)
    assert(strcmpi(opts.rel_chans, 'all'), 'Unrecognized input for rel_chans');
    opts.rel_chans = 1:numel(valid_data);
end

assert(~isempty(opts.rel_chans), 'No channels selected');
if ischar(opts.comps)
    assert(strcmpi(opts.comps, 'all'), 'Unrecognized input for comps');
    opts.comps = 1:size(valid_data{1}, 1);
end

% allow xcorr_pairs to be specified as an n x 2 matrix
if ~iscell(opts.xcorr_pairs)
    opts.xcorr_pairs = num2cell(opts.xcorr_pairs, 2);
end

Fs = 1 / data_s.options.winstep;  % window step = sample period of TFR data
sm_span_samp = Fs * opts.smooth_span;
opts.plot_raw = opts.plot_raw && sm_span_samp > 1; % if not smoothing, don't need to plot raw data

%-----------------------------------

% select channels and components
names = valid_names(opts.rel_chans);
data_raw = valid_data(opts.rel_chans);
data_raw = cellfun(@(pcd) pcd(opts.comps, :), data_raw, 'uni', false);

% smooth before diffing
data_sm = cellfun(@(pcd) smoothdata(pcd, 2, opts.smooth_method, sm_span_samp, 'includenan'), ...
    data_raw, 'uni', false);

% compute norm of change in pca data == change speed
pca_change_raw = cellfun(@(pcd) Fs * vecnorm(diff(pcd, 1, 2), opts.norm_type), data_raw, 'uni', false);
pca_change_raw = vertcat(pca_change_raw{:});

pca_change = cellfun(@(pcd) Fs * vecnorm(diff(pcd, 1, 2), opts.norm_type), data_sm, 'uni', false);
pca_change = vertcat(pca_change{:});

pca_time = data_s.time_grid;
pca_change_time = mean([pca_time(1:end-1); pca_time(2:end)]);

% interpolate if needed
if opts.interp_factor > 1
    time_interp = pca_change_time(1):1/(Fs*opts.interp_factor):pca_change_time(end);
    
    pca_change_raw_interp = interp1(pca_change_time, pca_change_raw.', time_interp, opts.interp_method).';
    pca_change = interp1(pca_change_time, pca_change.', time_interp, opts.interp_method).';
    pca_change_time = time_interp;
else
    pca_change_raw_interp = pca_change_raw;
end

% save, if requested
if opts.save
    data_s.(opts.name) = pca_change;
    data_s.([opts.name, '_time']) = pca_change_time;
    data_s.([opts.name, '_opts']) = opts;
    data_s.([opts.name, '_name']) = names;
    save(result_file, '-struct', 'data_s', '-v7.3');
end

% make plots, if requested
if opts.plot
    
    figure;
    for kC = 1:n_chans
        subplot(n_chans, 1, kC);
        
        if opts.plot_raw
            plot(pca_change_time, pca_change_raw_interp(kC, :), 'k');
            hold on;
        end
        
        plot(pca_change_time, pca_change(kC, :), 'b');
        axis tight;
        xlabel('Time (s)');
        ylabel('Euclidean change per second');
        title(sprintf('%s PCA-space change speed', names{kC}));
        
        if opts.plot_raw
            legend('Raw', sprintf('Smoothed (span = %d s)', opts.smooth_span));
        end        
    end
    
    if opts.savefigs
        savefig(fullfile(fileparts(result_file), 'pca_change.fig'));
    end
    
    n_xcorrs = length(opts.xcorr_pairs);
    if n_xcorrs > 0
        figure;
        
        dt = 1/Fs; % use non-interpolated version
        numlags = 20 / min(dt, 1);
        
        for kX = 1:n_xcorrs
            ax = subplot(n_xcorrs, 1, kX);
            
            pair = opts.xcorr_pairs{kX};
            crosscorr(pca_change_raw(pair(1), :), pca_change_raw(pair(2), :), 'NumLags', numlags);
            ax.XTickLabel = arrayfun(@(l) num2str(l * dt), ...
                str2double(ax.XTickLabel), 'uni', false);
            xlabel('Lag (s)');
            title(sprintf('Sample Cross Correlation - %s vs. %s PCA-space change', ...
                data_s.name{opts.rel_chans(pair(1))}, data_s.name{opts.rel_chans(pair(2))}));
        end
        
        if opts.savefigs
            savefig(fullfile(fileparts(result_file), 'pca_change_xcorr.fig'));
        end
    end
end

end