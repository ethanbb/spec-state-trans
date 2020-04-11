function res_data = mt_preprocess(data_or_filename, options)
% Preprocess the output of multitaper for further analyses.
% This includes smoothing/filtering in time and frequency axes and normalization.
%
% First input can either be:
%  * a mt_res filename, in which case results are loaded from and saved to the file
%    and the cell of processed spectrograms is returned.
%  * a MatFile object, which is treated the same as a filename after loading
%  * a struct, in which case results are takend from the field $in_name and saved to $name,
%    and the whole updated struct is returned.
%
% Valid normalization types:
%   'none'        - do nothing
%   'log'         - do nothing, output in dB
%   'log_meansub' - output in dB with 0 mean across time
%   'log_z'       - output in dB with 0 mean and 1 SD across time
%   'rank'        - output rank-order scores (values in range [0, 1])

% Load as MatFile object if a filename
if ischar(data_or_filename)
    res_data = matfile(data_or_filename);
else
    res_data = data_or_filename;
end

% Default options:
opts = struct(...
    'in_name',          'pxx',      ... field name of input
    'name',             'pxx_pp',   ... field name of output
    'freq_sm_type',     'med',      ... type of frequency smoothing - 'none', 'med', 'exp', or
                                    ... a valid method input to smoothdata
    'freq_sm_span',     40,         ... span of frequency smoothing, in Hz
    'time_sm_type',     'exp',      ... type of temporal smoothing - same options as freq_sm_type
    'time_sm_span',     60,         ... span of temporal smoothing, in seconds
    'frac_power',       false,      ... whether to divide by total power across frequencies at each time
    'norm_type',        'log_z'     ... type of normalization - see header comment
    );

if nargin < 2 || isempty(options)
    options = struct;
end

opts_in = fieldnames(options);
for kO = 1:length(opts_in)
    opts.(opts_in{kO}) = options.(opts_in{kO});
end

pxx = res_data.(opts.in_name);
n_chans = numel(pxx);

pxx_pp = cell(size(pxx));

for kC = 1:n_chans
    
    % Step 1: filter/smooth
    
    % ...over frequencies
    pxx_pp{kC} = smooth_pxx(pxx{kC}, opts.freq_sm_type, opts.freq_sm_span, 1);
    
    % ...over time
    mt_opts = res_data.options;
    time_span_samp = round(opts.time_sm_span / mt_opts.winstep);
    segs = res_data.seg_windows;
    pxx_pp{kC} = smooth_pxx(pxx_pp{kC}, opts.time_sm_type, time_span_samp, 2, segs);
    
    % Step 2: use fraction of total power, if requested
    if opts.frac_power
        pxx_pp{kC} = pxx_pp{kC} ./ sum(pxx_pp{kC});
    end
    
    % Step 3: normalize
    switch opts.norm_type
        case 'none'
            % do nothing
        case 'log'
            pxx_pp{kC} = 10*log10(pxx_pp{kC});
        case 'log_meansub'
            pxx_pp{kC} = 10*log10(pxx_pp{kC});
            pxx_pp{kC} = pxx_pp{kC} - nanmean(pxx_pp{kC}, 2);
        case 'log_z'
            pxx_pp{kC} = normalize(log(pxx_pp{kC}), 2);
        
        case 'rank'
            % only operate on non-nan values
            valid_t = find(~any(isnan(pxx_pp{kC})));
            ranks = linspace(0, 1, length(valid_t));
            [~, sortinds] = sort(pxx_pp{kC}(:, valid_t), 2);
            
            for kF = 1:size(pxx_pp{kC}, 1)
                pxx_pp{kC}(kF, valid_t(sortinds(kF, :))) = ranks;
            end
    end   
end

res_data.(opts.name) = pxx_pp;

% if we're not dealing with a struct, just return the preprocessed results
if ~isstruct(data_or_filename)
    res_data = res_data.(opts.name);
end

end

function pxx = smooth_pxx(pxx, type, span_samp, dim, segs)
% do smoothing of specified type

switch type
    case 'none'
        return;
    case 'med'
        pxx = medfilt1(pxx, span_samp, [], dim);
    case 'exp'
        % construct IIR filter with correct timescale
        decay_2t = normpdf(2.5) / normpdf(0);
        [b_exp, a_exp] = exp_filter(span_samp, decay_2t);
        
        if dim == 2
            assert(nargin >= 5, 'segs input is required for exp smoothing along time');
            pxx = filtfilt_segs(b_exp, a_exp, pxx, segs, 2);
        else
            valid_cols = ~any(isnan(pxx));
            pxx(:, valid_cols) = filtfilt(b_exp, a_exp, pxx(:, valid_cols));
        end
    otherwise
        % assume it is a valid input for smoothdata
        pxx = smoothdata(pxx, dim, type, span_samp, 'includenan');
end
end