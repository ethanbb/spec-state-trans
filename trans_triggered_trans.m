% Given transition times in each channel based on the max NMF component, plot histograms of
% the transition times of one channel relative to transitions in another, within a certain radius.

sr_dirs = prepSR;

days = {
    '2020-01-30'
    '2020-01-31'
    '2020-02-06'
    '2020-03-05'
    '2020-03-06'
    '2020-03-10'
    '2020-03-11'
};

% analysis params
Fw = 10; % window frequency
radius_s = 60; % how far back and forward to look from each transition point

for kD = 1:length(days)
    
    % load class vectors to find transition times
    nmf_info = load(fullfile(sr_dirs.results, days{kD}, 'nmf_res.mat'), 'nmf_classes', 'chan_names');
    % repetitions 1 & 2
    classes_A = nmf_info.nmf_classes{1};
    classes_B = nmf_info.nmf_classes{2}; 
    
    % find what channels should actually be used, based on CSD
    csd_chans_V1_s = load(fullfile(sr_dirs.results, days{kD}, 'csd_V1.mat'), 'chan_names');
    csd_chans_MC_s = load(fullfile(sr_dirs.results, days{kD}, 'csd_MC.mat'), 'chan_names');
    layers = {'L2/3', 'L4', 'L5'};
    b_chan = [ismember(layers, csd_chans_V1_s.chan_names), ...
        ismember(layers, csd_chans_MC_s.chan_names)];
    
    chan_names = nmf_info.chan_names(b_chan);
    n_chans = length(chan_names);
    classes_A = classes_A(b_chan);
    classes_B = classes_B(b_chan);
    
    % get transition times for each class.
    b_trans_A = cellfun(@(cls) [false; diff(cls) ~= 0], classes_A, 'uni', false);
    b_trans_B = cellfun(@(cls) [false; diff(cls) ~= 0], classes_B, 'uni', false);
    rec_len = length(b_trans_A{1}) / Fw;
    
    % will align classes_B transitions to classes_A transition times. The latter should exclude
    % times within the "radius" of start and end of the recording to avoid edge falloff.
    trans_times_A = cellfun(@(bt) find(bt(1+radius_s*Fw:end-radius_s*Fw)) / Fw + radius_s, b_trans_A, 'uni', false);
    trans_times_B = cellfun(@(bt) find(bt) / Fw, b_trans_B, 'uni', false);
    
    fh = figure('Position', [1, 41, 1536, 747]); % happens to be full-screen on my laptop
    tl = tiledlayout(n_chans, n_chans);
    title(tl, ['Aligned state transitions on ', days{kD}]);
    xlabel(tl, 'Relative transition times (s)');
    
    for chanA = 1:n_chans
        nameA = chan_names{chanA}; 
        tAs = trans_times_A{chanA};

        for chanB = 1:n_chans
            nameB = chan_names{chanB};
            tBs = trans_times_B{chanB};
            
            nexttile;
            rel_trans_c = arrayfun(@(tA) tBs(abs(tBs - tA) < radius_s) - tA, tAs, 'uni', false);
            rel_trans = cell2mat(rel_trans_c);
            histogram(rel_trans);
            
            % make channel labels on border
            if chanB == 1
                ylabel(['near ', nameA, ' trans.'], 'Interpreter', 'none', 'FontSize', 12, ...
                    'FontWeight', 'bold', 'Rotation', 0, 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'middle');
            end
            
            if chanA == 1
                title([nameB, ' transitions'], 'Interpreter', 'none', 'FontSize', 12);
            end
        end
    end
    
    savefig(fh, fullfile(sr_dirs.results, days{kD}, 'rel_trans.fig'));
end
