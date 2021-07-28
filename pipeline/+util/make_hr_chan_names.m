function hr_chan_names = make_hr_chan_names(chan_names, spacing_or_locs)
% transform cell or string array of "Sup"/"L4"/"Inf" channel names into a more human-readable form
% spacing_or_locs = distance between selected channels (e.g. Sup2 to Sup1) in um
% Or, can be a struct of relative channel locations in L4 (per probe)

chan_names = string(chan_names(:));
region_layer = split(chan_names, '_');

regions = unique(region_layer(:, 1), 'stable');
layers_byreg = cell(length(regions), 1);
hr_chan_names = layers_byreg;

for kR = 1:length(regions)
    layers_byreg{kR} = region_layer(region_layer(:, 1) == regions(kR), 2);
    hr_chan_names{kR} = layers_byreg{kR};
end

if isstruct(spacing_or_locs)
    spacing_or_locs = struct2cell(spacing_or_locs);
    assert(length(spacing_or_locs) == length(layers_byreg) && ...
        all(cellfun(@length, spacing_or_locs) == cellfun(@length, layers_byreg)), ...
        'Wrong number of regions or channels/region in location struct');
end

L4_depth = 650;

for kR = 1:length(layers_byreg)
    l4_ind = find(layers_byreg{kR} == "L4");
    assert(length(l4_ind) <= 1, 'There can only be one layer 4 per region');

    if iscell(spacing_or_locs)
        % if applicable, shift relative locations such that L4 has the right depth
        if ~isempty(l4_ind)
            spacing_or_locs{kR} = spacing_or_locs{kR} - spacing_or_locs{kR}(l4_ind) + L4_depth;
        end
        hr_chan_names{kR} = string(spacing_or_locs{kR}(:));

    else % we just have spacing between channels
        for kC = 1:length(layers_byreg{kR})
            if kC == l4_ind
                hr_chan_names{kR}(kC) = L4_depth;
            else
                supnum = sscanf(layers_byreg{kR}(kC), 'Sup%d');
                infnum = sscanf(layers_byreg{kR}(kC), 'Inf%d');
                
                if ~isempty(supnum)
                    hr_chan_names{kR}(kC) = L4_depth - spacing_or_locs * supnum;
                elseif ~isempty(infnum)
                    hr_chan_names{kR}(kC) = L4_depth + spacing_or_locs * infnum;
                else
                    error('Layer name not understood - not L4, Sup or Inf');
                end
            end
        end
    end
    
    % add L4 hint if we're in V1
    if contains(regions(kR), "V1")
        hr_chan_names{kR}(l4_ind) = hr_chan_names{kR}(l4_ind) + " (L4)";
    end
end

hr_chan_names = vertcat(hr_chan_names{:});

end
