function fh = pick_csd_channels(csd_dir, layers_um, layer_names)
% Given a directory containing CSD result csd.fig, display a simple GUI to 
% identify layer 4 in V1 and other layers relative to it. Defaults to
% L2/3 (or supergranular), L4, L5 (or infragranular) and L5B (deeper point in infragranular).
% Can override by providing layers_um and layer_names. Given depths are only relative to L4 position.

if ~exist('layers_um', 'var') || isempty(layers_um)
    layers_um = [270, 580, 955, 1200];
end

if ~exist('layer_names', 'var') || isempty(layer_names)
    layer_names = {'L2/3', 'L4', 'L5', 'L5B'};
end

[layers_um, layer_ind] = sort(layers_um);
layer_names = layer_names(layer_ind);
l4_ind = find(strcmp(layer_names, 'L4'));
assert(numel(l4_ind) == 1, 'Exactly one layer must be L4');

% shift so layer 4 = 0
layers_um = layers_um - layers_um(l4_ind);

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
isV1 = @(fig) contains(fig.Children(2).Title.String, 'V1');

csd_figs = openfig(fullfile(csd_dir, 'csd.fig'));
fig_is_V1 = arrayfun(isV1, csd_figs);
close(csd_figs(~fig_is_V1));
fh = figure(csd_figs(fig_is_V1));
fh.WindowStyle = 'docked';
ah = gca;
ah.OuterPosition(3) = 0.95;
hold on;

% Zoom out so that channels not included on the CSD can be marked
n_pad_chans = chan_axis(1) - 1;
total_chans = length(chan_axis) + 2 * n_pad_chans;
spacing = mean(diff(depth_axis)); % space b/w electrode contacts in cm
layers_chans = round(layers_um / 1000 / spacing);

yyaxis left;
ylim([depth_axis(1) - n_pad_chans * spacing, depth_axis(end) + n_pad_chans * spacing]);

yyaxis right;
ylim([1, total_chans]);


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
chan_inds_v1 = l4_chan + layers_chans;
chan_inds_v1(chan_inds_v1 < 1) = nan;
chan_inds_v1(chan_inds_v1 > total_chans) = nan;

% plot layers on the figure and choose which ones to keep
x_lims = ah.XLim;
cb_x_pos = 0.92;
cb_y_pos = @(cind) ah.Position(2) + ah.Position(4) * (ah.YLim(2) - cind) / diff(ah.YLim);

for kC = length(chan_inds_v1):-1:1
    ind = chan_inds_v1(kC);
    if isnan(ind)
        % skip this one
        cboxes(kC) = uicontrol('Style', 'checkbox', 'Value', 0, 'Visible', false);
        continue;
    end
    
    plot(x_lims, ind * [1, 1], 'r', 'LineWidth', 2);
    cb = uicontrol('Style', 'checkbox', 'Value', 1, ...
        'Units', 'normalized', 'String', sprintf('%s (#%d)', layer_names{kC}, ind));
    cb.Position(1) = cb_x_pos;
    cb.Position(2) = cb_y_pos(chan_inds_v1(kC)) - cb.Position(4) / 2;
    
    if kC < length(chan_inds_v1) && ~isnan(chan_inds_v1(kC + 1))
        % make sure they don't overlap
        last_top = cboxes(kC + 1).Position(2) + cb.Position(4);
        cb.Position(2) = max(cb.Position(2), last_top + (cb.Position(4) * 0.2));
    end
    
    cboxes(kC) = cb;
end


    function save_callback(~, ~)
        chans_to_keep = [cboxes.Value] == [cboxes.Max];
        layer_names = layer_names(chans_to_keep);
        chan_inds_v1 = chan_inds_v1(chans_to_keep);
        
        % save to V1 and MC files
        s_csd.chan_names = layer_names;
        s_csd.chans = chan_inds_v1;
        
        % now save to MC file as well
        %chan_inds_mc = chan_inds_v1 + total_chans;
        % not doing that anymore since I changed the order
        s_csd = matfile(fullfile(csd_dir, 'csd_M1.mat'), 'Writable', true);
        s_csd.chan_names = layer_names;
        s_csd.chans = chan_inds_v1;
        
        
        if isvalid(fh)
            close(fh);
        end
    end

save_btn = uicontrol('String', 'Save', 'Units', 'Normalized', 'Callback', @save_callback);
save_btn.Position(1) = 0.9;

disp('On the figure, check the channels that look correct and press "Save"');

end
