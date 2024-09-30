function [clIDs,clInds] = getIDbyDepth(root,dmin,dmax,clType)
%% Returns the cluster IDs of units between depths
%
% Inputs:
% root      = root structure
% dmin      = double; min depth in um (aka deepest point or nearest to probe tip)
% dmax      = double; max depth in um (aka shallowest point or farthest from probe tip)
% clType    = string; 'all','good','mua','noise' specifying which subgroup
%
% Outputs:
% clIDs = 1xN matrix of cluster IDs between dmin and dmax on the probe
%
% Created 9/4/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    dmin
    dmax
    clType = 'good'
end

if strcmp(clType,'all')
    clIDs  = root.info.cluster_id(root.info.depth > dmin & root.info.depth < dmax);
elseif strcmp(clType,'good')
    subIDs = root.info.cluster_id(root.goodind);
    clIDs  = subIDs(root.info.depth(root.goodind) > dmin & root.info.depth(root.goodind) < dmax);
elseif strcmp(clType,'mua')
    subIDs = root.info.cluster_id(root.muaind);
    clIDs  = subIDs(root.info.depth(root.muaind) > dmin & root.info.depth(root.muaind) < dmax);
elseif strcmp(clType,'noise')
    subIDs = root.info.cluster_id(root.noiseind);
    clIDs  = subIDs(root.info.depth(root.noiseind) > dmin & root.info.depth(root.noiseind) < dmax);
end

if ~isempty(clIDs)
    for i = 1:length(clIDs)
        clInds(i) = find(root.info.cluster_id == clIDs(i));
    end
else
    clInds = [];
end

end
