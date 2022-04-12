function [num_chans, spacing] = get_probe_model_info(model_name)
% Given the model name of a probe or other recording device (e.g. 'H3', 'E64-500-20-70'),
% return the number of electrodes and spacing in um.

if startsWith(model_name, 'H3')
    % Cambridge Neurotech 64-channel
    num_chans = 64;
    spacing = 20;

elseif startsWith(model_name, 'H4')
    % Cambridge Neurotech 32-channel
    num_chans = 32;
    spacing = 25;

elseif startsWith(model_name, {'E64-500-20-60', 'E64-500-20-70'})
    % NeuroNexus ECoG array
    num_chans = 64;
    spacing = 500;
else
    error('ProektLab:unknownProbe', 'Name %s does not match a known probe model', model_name);
end

