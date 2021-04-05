function hr_chan_names = make_hr_chan_names(chan_names, spacing)
% transform cell or string array of "Sup"/"L4"/"Inf" channel names
% into a more human-readable form
% spacing = distance between selected channels (e.g. Sup2 to Sup1)
% in um

chan_names = string(chan_names(:));
region_layer = split(chan_names, '_');

hr_chan_names = chan_names;

for kC = 1:length(chan_names)
    layer_name = region_layer(kC, 2);
    if layer_name == "L4"
        hr_chan_names{kC} = strrep(hr_chan_names{kC}, '_', ' ');
    else
        supnum = sscanf(layer_name, 'Sup%d');
        infnum = sscanf(layer_name, 'Inf%d');
        
        if ~isempty(supnum)
            hr_chan_names{kC} = ['+', num2str(spacing * supnum), '{\mu}m'];
        elseif ~isempty(infnum)
            hr_chan_names{kC} = ['-', num2str(spacing * infnum), '{\mu}m'];
        else
            error('Layer name not understood - not L4, Sup or Inf');
        end
    end
end

end
