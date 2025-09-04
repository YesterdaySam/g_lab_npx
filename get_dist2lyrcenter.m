function [dist2center] = get_dist2lyrcenter(root,centers)
%% Calculates distance of all units to layer center
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   center = Double or Mx1 array of centers for M layers
%
% Outputs:
%   dist2center = Nx1 array of distances to center for all units in root
%
% Created 8/11/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    centers = 1     % If vector, use input as bin centers, otherwise calculate from root variable
end

if length(centers) > 1
    tmpCtr = centers;
else
    % tmpCtr = root.lfpinfo.lfpDepth(root.uPSDMax(1,:));    % Use LFP max
    tmpCtr = root.lyrbounds(1,:) + (root.lyrbounds(2,:) - root.lyrbounds(1,:))/2;   %Use layer bound center
end

nShanks = numel(unique(root.info.shankID));

dist2center = zeros(height(root.info),1);

for i = 1:nShanks
    shUnits = root.info.shankID == i-1;
    dist2center(shUnits) = root.info.depth(shUnits) - tmpCtr(i);
end

end