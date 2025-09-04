function [root] = get_dist2border(root,borderDists)
%% Calculates distance of all units to layer center and assigns to root
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   borderDists = Mx1 array of distances for M shanks
%
% Outputs:
%   root = updated root with root.info.bdrDist field
%
% Created 8/11/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

nShanks = numel(unique(root.info.shankID));

bdrDist = zeros(height(root.info),1);

for i = 1:nShanks
    shUnits = root.info.shankID == i-1;
    bdrDist(shUnits) = borderDists(i);
end

root.info.bdrDist = bdrDist;

end