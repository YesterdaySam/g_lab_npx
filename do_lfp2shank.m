function [root] = do_lfp2shank(root,sh,depth)
%% Adds back extraneous LFP channels from lfp to root
%
% Inputs:
%   root = root object with uPSDMax and lfpinfo fields
%   sh   = 1-indexed shank ID (e.g. shank 4 = 4)
%   depth = 
% Outputs:
%   root = root object with 
%
% Created 6/3/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

newInd = find(root.lfpinfo.lfpDepth == depth & root.lfpinfo.lfpShank == sh-1);

if isempty(newInd)
    disp('No exact match found for specified depth, using nearest depth instead')
    minDist = min(abs(root.lfpinfo.lfpDepth - 421));
    newDepth = root.lfpinfo.lfpDepth(find(abs(root.lfpinfo.lfpDepth - 421) == minDist,1));
    newInd = find(root.lfpinfo.lfpDepth == newDepth & root.lfpinfo.lfpShank == sh-1);
end

root.uPSDMax(2,sh) = newInd;

end