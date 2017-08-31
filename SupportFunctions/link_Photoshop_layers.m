function [  ] = link_Photoshop_layers( layer_names )
% Robert F Cooper 10-17-2014
%   This function takes in a cell array of names for photoshop to link
%   together. It does this by making the first layer in the list active,
%   then adding all other layers to the existing selection.

    if length(layer_names) > 1
        setActiveLayer(layer_names{1});

        for i=2:length(layer_names)

            addToSelection( layer_names{i} )
            
        end

        linkSelectedLayers();
    end
end

