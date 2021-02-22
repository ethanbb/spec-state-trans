function plot_transitions(transition_table, color_var, opts, line_opts)
% Make a raster plot of transitions in each channel in a given table from get_state_transitions
% or assess_transitions_globality.
%
% Syntax: plot_transitions(transition_table, [color_var]='sync_score', [...opts...], [...Line opts...])
%
% color_var:  Should be a string; if this is a variable in the transition table,
%             colors each line to indicate the value according to the VarColormapFn.
%
% 
% Options:
%   - 'ColorRange' [[0, 1]]:     Range of values for the colorbar (use data range if empty)
%   - 'VarLabel' ['Sync score']: Label for the colorbar, if color_var is nonempty.
%   - 'AxesBgColor':             'Color' property of axes (default is 0.65 gray if using color for
%                                lines, to make yellows visible)
%   - 'VarColormapFn' [@jet]:    Function to use to generate the color map for globality
%
% Takes any "line" name-value arguments to change properties of the plotted lines, including
% 'Parent' which can be used to specify the axis.

arguments
    transition_table table
    color_var (1,:) char {varMustBeInTable(transition_table, color_var)} ...
        = empty_if_var_not_in_table(transition_table, 'sync_score')
    opts.ColorRange double = [0, 1]
    opts.VarLabel (1,:) char = 'Sync score'
    opts.AxesBgColor = ([1, 1, 1] - (~isempty(color_var)) * 0.35)
    opts.VarColormapFn function_handle = @jet
    line_opts.?matlab.graphics.chart.primitive.Line
end

use_color = ~isempty(color_var);
if use_color
    color_var_vals = transition_table.(color_var);
    [var_colors, cmap, crange] = vals2colormap(color_var_vals, opts.VarColormapFn, opts.ColorRange);
end

n_chans = length(unique(transition_table.chan));
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
    
    if use_color
        h.Color = var_colors(kT, :);
    end
end

xlabel('Time (s)');
yticks(1:n_chans);
yticklabels(chan_names);

% add colorbar
if use_color
    colormap(ax, cmap);
    caxis(crange);
    cb = colorbar;
    cb.Label.String = opts.VarLabel;
end

end

function answer = empty_if_var_not_in_table(tt, var)
if ismember(var, tt.Properties.VariableNames)
    answer = var;
else
    answer = '';
end
end

% validation function
function varMustBeInTable(transition_table, var)
if ~isempty(var) && isempty(empty_if_var_not_in_table(transition_table, var))
    error('Cannot color based on %s - info not present', var);
end
end
