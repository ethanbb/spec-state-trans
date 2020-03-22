function res = multitaper_analysis(data_s, options)

% Do multitaper analysis of multiple LFP channels in the same way as
% Hudson et al., 2014
%
% data_s = struct from loading the raw dataset
% options = struct with overrides of the following defaults:
%
% If there are artifacts specified, the outputs pxx and phase (if
% requested) will have NaNs at timepoints that would be affected by the
% artifact.

%---------- organize data ----------%

Fs = data_s.finalSampR;
lfp_organized = organize_lfp(data_s);

%---------- parse arguments --------%
opts = struct(...
   'chans',      [16, 48],       ... channel numbers to analyze (default = middle channel from each region)
   'chan_names', {{'VC', 'MC'}}, ...
   'artifacts',  zeros(0, 2),    ... k x 2 list or k-length cell of [start, end] of artifacts, in seconds
   'n_tapers',   17,             ... number of tapers
   'window',     60,             ... window length (sec)
   'pad',        [],             ... pad length (samples, default = next power of 2)
   'padbase',    [],             ... pad windows as if they were this long, in seconds
   'winstep',    1,              ... window timestep (sec)
   'max_freq',   300,            ... max freq of interest (min freq is determined by n_tapers & window)
   'inc_phase',  false,          ... whether to include the phase
   'save',       true,           ... whether to save the results struct before returning
   'savedir',    pwd,            ... directory in which to save
   'filename',   'mt_res.mat'    ... filename to save
);

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

% derive pad from padbase if necessary
if ~isempty(opts.padbase)
    assert(isempty(opts.pad), 'Cannot specify both pad and padbase');
    opts.pad = 2.^max(nextpow2(opts.padbase * Fs), 14);
elseif isempty(opts.pad) % base padding on window
    opts.pad = 2.^max(nextpow2(opts.window * Fs),14);
end

% infer clean segments from artifacts
len_secs = size(lfp_organized, 2) / Fs;
if iscell(opts.artifacts)
    opts.artifacts = cellfun(@(a) a(:).', opts.artifacts, 'uni', false);
    opts.artifacts = vertcat(opts.artifacts{:});
end
opts.artifacts = sortrows(opts.artifacts);

opts.clean_segs = [[0; opts.artifacts(:, 2)], [opts.artifacts(:, 1); len_secs]];
clean_len = diff(opts.clean_segs, 1, 2);

% validate and remove empty segments
assert(all(clean_len >= 0), 'Overlapping or invalid artifacts');
opts.clean_segs(clean_len == 0, :) = [];

% compute minimum frequency from n_tapers and window size
thfp = (opts.n_tapers + 1) / 2;
opts.min_freq = thfp / (opts.pad / Fs);

n_chans = length(opts.chans);
n_segs = size(opts.clean_segs, 1);

% compute time grid for whole input (see lines 83-91 of swTFspecAnalog.m)
n_steps = round((len_secs - opts.window) / opts.winstep);
if len_secs == opts.window
    n_steps = 1;
end
res.time_grid = opts.window/2 + ((1:n_steps)-1)*opts.winstep;
win_mins = res.time_grid - opts.window/2;
win_maxes = res.time_grid + opts.window/2;

res.options = opts;
res.name = opts.chan_names(:);
res.pxx = cell(n_chans, 1);
if opts.inc_phase
    res.phase = cell(n_chans, 1);
end

% in case there are multiple segments, save what inputs and what outputs
% correspond to each one
res.seg_samples = cell(1, n_segs);
res.seg_windows = cell(1, n_segs);


for kS = 1:n_segs
    
    %------ identify which windows to process ------%
    
    res.seg_windows{kS} = find(win_mins >= opts.clean_segs(kS, 1) & ...
                               win_maxes <= opts.clean_segs(kS, 2));
    
    res.seg_samples{kS} = win_mins(res.seg_windows{kS}(1))*Fs + 1 : ...
                          win_maxes(res.seg_windows{kS}(end))*Fs;
    
    for kC = 1:n_chans
               
        %------ grab the segment of data -------%
        
        % pad with zeros to force it to process all the windows
        extra_zeros = zeros(1, opts.winstep * Fs);
        
        lfp_seg = [lfp_organized(opts.chans(kC), res.seg_samples{kS}), extra_zeros];
                
        %------ do multitaper analysis ---------%
        
        mt_res = swTFspecAnalog(lfp_seg, Fs, opts.n_tapers, ...
            [opts.min_freq, opts.max_freq], opts.window * Fs, ...
            opts.winstep * Fs, [], opts.pad, [], [], [], opts.inc_phase);
        
        if kS == 1
            if kC == 1
                % on the very first call, make frequency grid
                res.freq_grid = mt_res.freq_grid;
                n_freqs = length(res.freq_grid);
            end
            
            % initialize all to nans, then fill in the steps that get
            % computed
            res.pxx{kC} = nan(n_freqs, n_steps);
            if opts.inc_phase
                res.phase{kC} = nan(n_freqs, n_steps);
            end
        end
        
        res.pxx{kC}(:, res.seg_windows{kS}) = squeeze(mt_res.tfse).';
        if opts.inc_phase
            res.phase{kC}(:, res.seg_windows{kS}) = squeeze(mt_res.phase).';
        end
    end
end

if opts.save
    if ~exist(opts.savedir, 'dir')
        % attempt to create it
        mkdir(opts.savedir);

        if ~exist(opts.savedir, 'dir')
            warning('Directory %s could not be created; saving in current dir instead', ...
                opts.savedir);
            opts.savedir = pwd;
        end
    end

    save(fullfile(opts.savedir, opts.filename), '-struct', 'res', '-v7.3');
end

end