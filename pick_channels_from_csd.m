% For all recording days, use CSD analysis and diagram of rat cortical layers
% to pick channels to analyze. Aiming for 1 superficial (layer 2/3), 1 layer 4, and 1 layer 5.

% Get recording days

prepSR;

days = {
    '2020-01-30'
    '2020-01-31'
    '2020-02-06'
    '2020-03-05'
    '2020-03-06'
    '2020-03-10'
    '2020-03-11'
    };
n_days = length(days);

csd_dirs = fullfile(results_dir, days);

for kD = 1:n_days
    uiwait(pick_channels(csd_dirs{kD}));
end

function fh = pick_channels(csd_dir)

% positions of layers based on diagram w/ surface = 0 (2/3, 4, 5)
layers_um = [270, 580, 955];

% shift so layer 4 = 0
layers_um = layers_um - layers_um(2);

% convert to relative number of channels
spacing = 25; % space b/w electrode contacts in um
layers_chans = round(layers_um / spacing);

% Now locate layer 4 on the V1 CSD (by just finding the peak)
s_csd = matfile(fullfile(csd_dir, 'csd_V1.mat'), 'Writable', true);
time_axis = s_csd.time_axis;
depth_axis = s_csd.depth_axis;
chan_axis = s_csd.chan_axis;
csd = s_csd.csd;

b_after_stim = time_axis >= 0 & time_axis <= 100;
time_after_stim = time_axis(b_after_stim);
csd_after_stim = csd(:, b_after_stim);
[max_across_chans, maxind_across_chans] = max(csd_after_stim);
[~, maxind_across_time] = max(max_across_chans);
max_time = time_after_stim(maxind_across_time);
l4_chan = s_csd.chan_axis(1, maxind_across_chans(maxind_across_time));

% Open the figure to allow verification/manual correction if necessary.
csd_figs = openfig(fullfile(csd_dir, 'csd.fig'));
close(csd_figs(2));
fh = figure(csd_figs(1));
fh.WindowStyle = 'docked';
ah = gca;
ah.OuterPosition(3) = 0.95;
hold on;

% Zoom out so that channels not included on the CSD can be marked
n_extra_start = chan_axis(1) - 1;
n_extra_end = 32 - chan_axis(end);
spacing = mean(diff(depth_axis));

yyaxis left;
ylim([depth_axis(1) - n_extra_start * spacing, depth_axis(end) + n_extra_end * spacing]);

yyaxis right;
ylim([1, 32]);


plot(max_time, l4_chan, 'r.', 'MarkerSize', 15);

while true
    new_chan = input('Does the marked spot look like L4? If so, press enter; if not, enter which channel it should be: ');
    if isempty(new_chan) || (isscalar(new_chan) && ismember(new_chan, s_csd.chan_axis))
        break;
    else
        disp('Invalid input - either enter nothing or the channel layer 4 should be.');
    end
end

if ~isempty(new_chan)
    l4_chan = new_chan;
end

% save best possible chans for V1 to matfile
chan_names = {'L2/3', 'L4', 'L5'};
chan_inds_v1 = max(1, min(32, l4_chan + layers_chans));

% plot layers on the figure and choose which ones to keep
x_lims = ah.XLim;
cb_x_pos = 0.92;
cb_y_pos = @(cind) ah.Position(2) + ah.Position(4) * (ah.YLim(2) - cind) / diff(ah.YLim);

for kC = length(chan_inds_v1):-1:1
    ind = chan_inds_v1(kC);
    plot(x_lims, ind * [1, 1], 'r', 'LineWidth', 2);
    cb = uicontrol('Style', 'checkbox', 'Value', 1, ...
        'Units', 'normalized', 'String', sprintf('%s (#%d)', chan_names{kC}, ind));
    cb.Position(1) = cb_x_pos;
    cb.Position(2) = cb_y_pos(chan_inds_v1(kC)) - cb.Position(4) / 2;
    
    if kC < length(chan_inds_v1)
        % make sure they don't overlap
        last_top = cboxes(kC + 1).Position(2) + cb.Position(4);
        cb.Position(2) = max(cb.Position(2), last_top + (cb.Position(4) * 0.2));
    end
    
    cboxes(kC) = cb;
end


    function save_callback(~, ~)
        chans_to_keep = [cboxes.Value] == [cboxes.Max];
        chan_names = chan_names(chans_to_keep);
        chan_inds_v1 = chan_inds_v1(chans_to_keep);
        
        % save to V1 and MC files
        s_csd.chan_names = chan_names;
        s_csd.chans = chan_inds_v1;

        % now save to MC file as well
        chan_inds_mc = chan_inds_v1 + 32;
        s_csd = matfile(fullfile(csd_dir, 'csd_MC.mat'), 'Writable', true);
        s_csd.chan_names = chan_names;
        s_csd.chans = chan_inds_mc;

        
        if isvalid(fh)
            close(fh);
        end
    end

save_btn = uicontrol('String', 'Save', 'Units', 'Normalized', 'Callback', @save_callback);
save_btn.Position(1) = 0.9;

disp('On the figure, check the channels that look correct and press "Save"');

end