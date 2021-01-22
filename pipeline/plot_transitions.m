function plot_transitions(transition_table, color_globality, opts, line_opts)
% Make a raster plot of transitions in each channel in a given table from get_state_transitions
% or assess_transitions_globality.
%
% Syntax: plot_transitions(transition_table, [color_globality]=true, [...opts...], [...Line opts...])
%
% color_globality: if true (the default), colors each line to indicate the globality of the
%                  transition (see assess_transition_globality.m).
% 
% Options:
%   - 'AxesBgColor': 'Color' property of axes (default is 0.65 gray if color_globality is true,
%                     to make yellows visible)
%   - 'GlobalityColormapFn': Function to use to generate the color map for globality
%
% Takes any "line" name-value arguments to change properties of the plotted lines, including
% 'Parent' which can be used to specify the axis.

arguments
    transition_table table
    color_globality (1,1) logical {mustBeFalseIfNoGlobality(color_globality, transition_table)} ...
        = has_globality(transition_table)
    opts.AxesBgColor = ([1, 1, 1] - color_globality * 0.35)
    opts.GlobalityColormapFn function_handle = @jet
    line_opts.?matlab.graphics.chart.primitive.Line
end

n_chans = length(unique(transition_table.chan));
globality_colors = opts.GlobalityColormapFn(n_chans);

chan_names = arrayfun(@(c) unique(transition_table.chan_name(transition_table.chan == c)), 1:n_chans);

if isfield(line_opts, 'Parent')
    axes(line_opts.Parent);
end
ax = gca;
ax.Color = opts.AxesBgColor;
ax.TickLabelInterpreter = 'none';
ax.YDir = 'reverse';
hold on;

line_opts_cell = namedargs2cell(line_opts);
for kT = 1:height(transition_table)
    chan = transition_table.chan(kT);
    time = transition_table.time(kT);
    h = plot([time, time], [chan-0.5, chan+0.5], 'k', line_opts_cell{:});
    
    if color_globality
        h.Color = globality_colors(transition_table.globality(kT), :);
    end
end

xlabel('Time (s)');
yticks(1:n_chans);
yticklabels(chan_names);

% add colorbar
if color_globality
    colormap(globality_colors);
    caxis([0.5, n_chans + 0.5]);
    cb = colorbar('Ticks', 1:n_chans);
    cb.Label.String = 'Globality (channels)';
end

end

function answer = has_globality(tt)
answer = ismember('globality', tt.Properties.VariableNames);
end

% validation function
function mustBeFalseIfNoGlobality(color_globality, transition_table)
if color_globality && ~has_globality(transition_table)
    error('Cannot color based on globality - info not present');
end
end
