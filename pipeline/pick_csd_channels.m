function fh = pick_csd_channels(csd_dir, layers_um, layer_names, ...
        main_probe, aligned_probes, dock_figs)
% Given a directory containing CSD result csd.fig, display a simple GUI to 
% identify layer 4 in V1 and other layers relative to it. Defaults to
% L2/3 (or supergranular), L4, L5 (or infragranular) and L5B (deeper point in infragranular).
% Can override by providing layers_um and layer_names. Given depths are only relative to L4 position.
% main_probe is the probe name/region to plot CSD of and select channels for (e.g. 'V1')
% aligned_probes is a cell of other probe names/regions which are aligned and into whose matfiles the
% same channel numbers should also be saved (e.g. {'M1'})
%
% In a script, to wait until the selection process is done before proceeding, call with uiwait,
% as in: uiwait(pick_csd_channels(...));

if ~exist('layers_um', 'var') || isempty(layers_um)
    layers_um = [270, 580, 955, 1200];
end

if ~exist('layer_names', 'var') || isempty(layer_names)
    layer_names = {'L2/3', 'L4', 'L5', 'L5B'};
end

if exist('aligned_probes', 'var')
    assert(iscell(aligned_probes), 'Aligned probes must be passed as a cell');
else
    aligned_probes = {};
end

if ~exist('dock_figs', 'var') || isempty(dock_figs)
    dock_figs = false;
end

[layers_um, layer_ind] = sort(layers_um);
layer_names = layer_names(layer_ind);
l4_ind = find(strcmp(layer_names, 'L4'));
assert(numel(l4_ind) == 1, 'Exactly one layer must be L4');

% shift so layer 4 = 0
layers_um = layers_um - layers_um(l4_ind);

% Now locate layer 4 on the V1 CSD (by just finding the peak)
s_csd = matfile(fullfile(csd_dir, sprintf('csd_%s.mat', main_probe)), 'Writable', true);
time_axis = s_csd.time_axis;
depth_axis = s_csd.depth_axis;
chan_axis = s_csd.chan_axis;
csd = s_csd.csd;
bad_chans = chan_axis(any(isnan(csd), 2));

b_after_stim = time_axis >= 0 & time_axis <= 100;
time_after_stim = time_axis(b_after_stim);
csd_after_stim = csd(:, b_after_stim);
[max_across_chans, maxind_across_chans] = max(csd_after_stim);
[~, maxind_across_time] = max(max_across_chans);
max_time = time_after_stim(maxind_across_time);
l4_chan = s_csd.chan_axis(1, maxind_across_chans(maxind_across_time));

% Open the figure to allow verification/manual correction if necessary.
isV1 = @(fig) contains(fig.Children(2).Title.String, main_probe);

csd_figs = openfig(fullfile(csd_dir, 'csd.fig'));
fig_is_V1 = arrayfun(isV1, csd_figs);
close(csd_figs(~fig_is_V1));
fh = figure(csd_figs(fig_is_V1));
if dock_figs
    fh.WindowStyle = 'docked';
end
ah = gca;
ah.OuterPosition(3) = 0.95;
x_lims = ah.XLim;
hold on;

for c = bad_chans(:)'
    bad_h = plot(ah, x_lims, c * [1, 1], 'r-', 'LineWidth', 1.5);
end

grid on;
ah.Layer = 'top'; % make gridlines visible

% Zoom out so that channels not included on the CSD can be marked
n_pad_chans = chan_axis(1) - 1;
total_chans = length(chan_axis) + 2 * n_pad_chans;
spacing = mean(diff(depth_axis)); % space b/w electrode contacts in cm
layers_chans = round(layers_um / 1000 / spacing);

ylim([1, total_chans]);

plot(max_time, l4_chan, 'r.', 'MarkerSize', 15);
if exist('bad_h', 'var')
    legend(bad_h, 'Bad channels', 'Location', 'northeast');
end

while true
    fprintf('Does the marked spot (channel %d) look like L4?\n', l4_chan);
    new_chan = input('If so, press enter; if not, enter which channel it should be: ');
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
chan_inds_v1 = l4_chan + layers_chans;
chan_inds_v1(chan_inds_v1 < 1) = nan;
chan_inds_v1(chan_inds_v1 > total_chans) = nan;

% deal with if any are bad chans
selected_bad = find(ismember(chan_inds_v1, bad_chans));
for kB = 1:length(selected_bad)
    bad_ind = selected_bad(kB);
    bad_chan = chan_inds_v1(bad_ind);
    
    % try to substitute with one lower or higher
    bad_chans_aug = [0; bad_chans(:); total_chans + 1];
    
    if ~ismember(bad_chan + 1, bad_chans_aug)
        chan_inds_v1(bad_ind) = bad_chan + 1;
    elseif ~ismember(bad_chan - 1, bad_chans_aug)
        chan_inds_v1(bad_ind) = bad_chan - 1;
    else
        % just skip it
        chan_inds_v1(bad_ind) = nan;
    end
end

% plot layers on the figure and choose which ones to keep
cb_x_pos = 0.9;
cb_width = 0.09;
cb_y_pos = @(cind) ah.Position(2) + ah.Position(4) * (ah.YLim(2) - cind) / diff(ah.YLim);

for kC = length(chan_inds_v1):-1:1
    ind = chan_inds_v1(kC);
    if isnan(ind)
        % skip this one
        cboxes(kC) = uicontrol(fh, 'Style', 'checkbox', 'Value', 0, 'Visible', false);
        continue;
    end
    
    sel_h = plot(ah, x_lims, ind * [1, 1], 'g-', 'LineWidth', 1.5);
    cb = uicontrol(fh, 'Style', 'checkbox', 'Value', 1, ...
        'Units', 'normalized', 'String', sprintf('%s (#%d)', layer_names{kC}, ind));
    cb.Position(1) = cb_x_pos;
    cb.Position(3) = cb_width;
    cb.Position(2) = cb_y_pos(chan_inds_v1(kC)) - cb.Position(4) / 2;
    
    if kC < length(chan_inds_v1) && ~isnan(chan_inds_v1(kC + 1))
        % make sure they don't overlap
        last_top = cboxes(kC + 1).Position(2) + cb.Position(4);
        cb.Position(2) = max(cb.Position(2), last_top + (cb.Position(4) * 0.2));
    end
    
    cboxes(kC) = cb;
end

if exist('bad_h', 'var')
    legend([bad_h, sel_h], {'Bad channels', 'Selected'}, 'Location', 'northeast');
else
    legend(sel_h, {'Selected'}, 'Location', 'northeast');
end

    function save_callback(~, ~)
        chans_to_keep = [cboxes.Value] == [cboxes.Max];
        layer_names = layer_names(chans_to_keep);
        chan_inds_v1 = chan_inds_v1(chans_to_keep);
        
        s_csd.chan_names = layer_names;
        s_csd.chans = chan_inds_v1;
        
        % now save to aligned probe files as well
        for kP = 1:numel(aligned_probes)
            probe = aligned_probes{kP};
            s_csd = matfile(fullfile(csd_dir, sprintf('csd_%s.mat', probe)), 'Writable', true);
            s_csd.chan_names = layer_names;
            s_csd.chans = chan_inds_v1;
        end
        
        if isvalid(fh)
            close(fh);
        end
    end

save_btn = uicontrol(fh, 'String', 'Save', 'Units', 'Normalized', 'Callback', @save_callback);
save_btn.Position(1) = 0.9;

disp('On the figure, check the boxes for the channels to keep for analysis and press "Save"');
uiwait(fh);

end
