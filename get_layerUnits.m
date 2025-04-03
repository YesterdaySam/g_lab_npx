function [root] = get_layerUnits(root)
%% Assigns units based on peak in ripple band power
%
% Inputs:
% root = root object. Must have lyrBounds field
%
% Outputs:
% root = updated root with root.lyrbounds and root.info.lyrID where:
%   0 = below layer
%   1 = in layer
%   2 = above layer
%
% Created 3/11/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

%% Assign units to layers by depth
nUnits = height(root.info);
for i = 1:nUnits
    tmpSh = root.info.shankID(i)+1;
    if root.info.depth(i) < root.lyrbounds(1,tmpSh)
        root.info.lyrID(i) = 0;
    elseif root.info.depth(i) <= root.lyrbounds(2,tmpSh)
        root.info.lyrID(i) = 1;
    else
        root.info.lyrID(i) = 2;
    end
end

end