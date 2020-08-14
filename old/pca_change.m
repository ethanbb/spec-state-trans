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
    'smooth_method', 'gaussian',     ... smoothing method (3rd argument of smoothdata).
                                     ... can also be 'exp' for an exponential (Poisson) window.
                                     ... see: https://www.desmos.com/calculator/hkwshtnbfn
    'diff_factor',   0.75,           ... step by 'smooth_span' times this factor when calculating finite difference
                                     ... see https://www.desmos.com/calculator/hhnuua3vau - factor 'a'
    'diff_step',     [],             ... overrides diff_factor to always use a certain diff step (in seconds)
    'norm_type',     2,              ... second input to vecnorm (default is Euclidian distance)
                                     ... or 'none' to output velocity (components in dimension 3)
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

Fw = 1 / data_s.options.winstep;  % window step = sample period of TFR data
sm_span_samp = Fw * opts.smooth_span;
opts.plot_raw = opts.plot_raw && sm_span_samp > 1; % if not smoothing, don't need to plot raw data

%-----------------------------------

% select channels and components
names = valid_names(opts.rel_chans);
data_raw = valid_data(opts.rel_chans);
data_raw = cellfun(@(pcd) pcd(opts.comps, :), data_raw, 'uni', false);

% smooth before diffing
if strcmpi(opts.smooth_method, 'exp')

    % manually implement exponential smoothing
    % match width for same decay @ edge of Gaussian w/ 2.5 SDs (default for smoothdata)
    decay_2t = normpdf(2.5) / normpdf(0);
    [b_exp, a_exp] = exp_filter(sm_span_samp, decay_2t);
    
    data_sm = filtfilt_segs(b_exp, a_exp, data_raw, data_s.seg_windows, 2);
else
    data_sm = cellfun(@(pcd) smoothdata(pcd, 2, opts.smooth_method, sm_span_samp, 'includenan'), ...
        data_raw, 'uni', false);
end

if ~isempty(opts.diff_step)
    diff_kernel_samps = opts.diff_step * Fw;
else
    % implement diff with step as in https://www.desmos.com/calculator/hhnuua3vau
    diff_kernel_samps = round(sm_span_samp * opts.diff_factor);
end
diff_kernel_samps = diff_kernel_samps + 1; % actual length of kernel

diff_kernel = [1, zeros(1, diff_kernel_samps - 2), -1];
mydiff = @(pcd) convn(pcd, diff_kernel, 'valid');

% function to diff and either take a norm or leave the velocity
scale_factor = Fw / (diff_kernel_samps - 1);
if strcmpi(opts.norm_type, 'none')
    changefn = @(pcd) scale_factor * permute(mydiff(pcd), [3, 2, 1]);
else
    changefn = @(pcd) scale_factor * vecnorm(mydiff(pcd), opts.norm_type, 1);
end

% compute norm of change in pca data == change speed
pca_change_raw = cellfun(changefn, data_raw, 'uni', false);
pca_change_raw = vertcat(pca_change_raw{:});

pca_change = cellfun(changefn, data_sm, 'uni', false);
pca_change = vertcat(pca_change{:});

pca_time = data_s.time_grid(:).';
pca_change_time = convn(pca_time, ones(1, diff_kernel_samps) / diff_kernel_samps, 'valid');

% interpolate if needed
if opts.interp_factor > 1
    pca_change_time = pca_change_time(1):1/(Fw*opts.interp_factor):pca_change_time(end);
    
    pca_change_raw_interp = interp1(pca_change_time, permute(pca_change_raw, [2, 1, 3]), ...
        pca_change_time, opts.interp_method);
    pca_change_raw_interp = permute(pca_change_raw_interp, [2, 1, 3]);
    
    pca_change = interp1(pca_change_time, permute(pca_change, [2, 1, 3]), ...
        pca_change_time, opts.interp_method);
    pca_change = permute(pca_change, [2, 1, 3]);
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
            plot(pca_change_time, pca_change_raw_interp(kC, :, 1), 'k');
            hold on;
        end
        
        plot(pca_change_time, pca_change(kC, :, 1), 'b');
        if size(pca_change, 3) > 1
            warning('Plotting only first component of change velocity');
        end
        
        axis tight;
        xlabel('Time (s)');
        ylabel('Change per second');
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
        
        dt = 1/Fw; % use non-interpolated version
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
