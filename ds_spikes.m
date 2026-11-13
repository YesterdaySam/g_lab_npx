function [root2, propDiff] = ds_spikes(root,dsProp)
%% 
% Inputs:
%   root = root object
%   dsProp = proportion of spikes to save
%
% Outputs:
%   root2 = updated root with fewer root.good spikes
%   propDiff = gut check - number of remaining root.good spikes / original
%   spikes, should be as close to 0 as possible
%
% Updated 10/30/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    dsProp
end

root2 = root;
nSpkOrig = nan(size(root.good));
rmArray = false(size(root.cl));
rootind = 1:length(root.cl);

for i = 1:length(root.good)
    tmpInds = root.cl == root.good(i);
    nSpkOrig(i) = sum(tmpInds);
    rmNum = round(nSpkOrig(i)*(1-dsProp));
    rmInds = datasample(rootind(tmpInds),rmNum,'Replace',false);
    rmArray(rmInds) = true;
end

root2.cl(rmArray) = [];
root2.ts(rmArray) = [];
root2.tsb(rmArray) = [];

propDiff = sum(rmArray) / sum(nSpkOrig) - (1-dsProp);  % Gut check
end