function [root] = assignRegion(root, regions)
%% Quick and dirty assignment of units by shank
%
% Inputs:
% root = root object.
%
% Outputs:
% root = root object with new subfield of root.info.region
%
% Created 1/24/2025 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

if root.prbType == 'NPX2.0'
    for i = 0:3
        tmpUnits = root.info.shankID == i;
        root.info.region(tmpUnits) = regions(i+1);
    end
else
    disp(['No assignment for 1.0 probe'])
end

end