function output_txt = scatter_update(obj,event_obj)
% Display data cursor position in a data tip
% obj          Currently not used
% event_obj    Handle to event object
% output_txt   Data tip text, returned as a character vector or a cell array of character vectors

pos = event_obj.Position;
series = event_obj.Target;
xinds = find(series.XData == pos(1));
yinds = find(series.YData == pos(2));
ind = intersect(xinds, yinds);


%********* Define the content of the data tip here *********%

% Display the x and y values:
output_txt = {['X',formatValue(pos(1),event_obj)],...
              ['Y',formatValue(pos(2),event_obj)]};


if isprop(series, 'SizeData') && length(series.SizeData) == length(series.XData)
    output_txt{end+1} = ['Size', formatValue(series.SizeData(ind(1)), event_obj)];
end

if isprop(series, 'CData') && length(series.CData) == length(series.XData)
    output_txt{end+1} = ['Color', formatValue(series.CData(ind(1)), event_obj)];
end
    
%***********************************************************%


% If there is a z value, display it:
if length(pos) > 2
    output_txt{end+1} = ['Z',formatValue(pos(3),event_obj)];
end

%***********************************************************%

function formattedValue = formatValue(value,event_obj)
% If you do not want TeX formatting in the data tip, uncomment the line below.
% event_obj.Interpreter = 'none';
if strcmpi(event_obj.Interpreter,'tex')
    valueFormat = ' \color[rgb]{0 0.6 1}\bf';
    removeValueFormat = '\color[rgb]{.25 .25 .25}\rm';
else
    valueFormat = ': ';
    removeValueFormat = '';
end
formattedValue = [valueFormat num2str(value,4) removeValueFormat];
