% For all recording days, use CSD analysis and diagram of rat cortical layers
% to pick channels to analyze. Aiming for 1 superficial (layer 2/3), 1 layer 4, and 1 layer 5.

% Get recording days

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

for kD = 1:n_days
    pick_channels(days{kD});
end

function pick_channels(day)

prepSR;

% positions of layers based on diagram w/ surface = 0 (2/3, 4, 5)
layers_um = [270, 580, 955];

% shift so layer 4 = 0
layers_um = layers_um - layers_um(2);

% convert to relative number of channels
spacing = 25; % space b/w electrode contacts in um
layers_chans = round(layers_um / spacing);

% Now locate layer 4 on the V1 CSD (by just finding the peak)
s_csd = matfile(fullfile(results_dir, day, 'csd_V1.mat'), 'Writable', true);
time_axis = s_csd.time_axis;
csd = s_csd.csd;

b_after_stim = time_axis >= 0 & time_axis <= 100;
time_after_stim = time_axis(b_after_stim);
csd_after_stim = csd(:, b_after_stim);
[max_across_chans, maxind_across_chans] = max(csd_after_stim);
[~, maxind_across_time] = max(max_across_chans);
max_time = time_after_stim(maxind_across_time);
l4_chan = s_csd.chan_axis(1, maxind_across_chans(maxind_across_time));

% Open the figure to allow verification/manual correction if necessary.
csd_figs = openfig(fullfile(results_dir, day, 'csd.fig'));
close(csd_figs(2));
figure(csd_figs(1));
hold on;
yyaxis right;
plot(max_time, l4_chan, 'r.', 'MarkerSize', 15);

while true
    new_chan = input('Does the marked spot look like L4? If so, press enter; if not, enter which channel it should be: ');
    if isempty(new_chan) || (isscalar(new_chan) && ismember(new_chan, s_csd.chan_axis))
        break;
    else
        disp('Invalid input - either enter nothing or the channel layer 4 should be.');
    end
end

if isvalid(csd_figs(1))
    close(csd_figs(1));
end

if ~isempty(new_chan)
    l4_chan = new_chan;
end

% save best possible chans for V1 to matfile
chan_names = {'L2/3', 'L4', 'L5'};
chan_inds_v1 = max(1, min(32, l4_chan + layers_chans));

s_csd.chan_names = chan_names;
s_csd.chans = chan_inds_v1;

% now save to MC file as well
chan_inds_mc = chan_inds_v1 + 32;
s_csd = matfile(fullfile(results_dir, day, 'csd_MC.mat'), 'Writable', true);
s_csd.chan_names = chan_names;
s_csd.chans = chan_inds_mc;

end