function [root] = burst2root(root,sess,thresh)
%% Finds bursts under thresh for each unit in root, adds burst onset time indices, cluster ID, and length
%
% Bursts >2 spikes in <6ms interleaved Simonnet et al 2019 (Juxtacellular patch) https://www.jneurosci.org/content/39/19/3651
% Bursts >2 spikes in <9ms Royer et al 2016 https://pmc.ncbi.nlm.nih.gov/articles/PMC4919905/
%
% Inputs:
%   root = root struct
%   sess = session struct
%   thresh = theshold in seconds for burst inter-spike-interval maximum
%
% Outputs:
%   root = updated root with burst_tsb, burst_cl and burst_len fields
%       burst_tsb = sorted index of behavioral timestamps for sess
%       burst_cl = cluster ID of each burst
%       burst_len = length of spikes in each burst
%
% Created 12/5/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root struct
    sess struct
    thresh = 0.008
end

root.burst_tsb = [];
root.burst_cl = [];
root.burst_len = [];

for i = 1:height(root.info)
    cc = root.info.cluster_id(i);
    [bStts,~,bLens] = get_bursts(root,sess,cc,thresh);

    root.burst_tsb = [root.burst_tsb; bStts];
    root.burst_cl = [root.burst_cl; cc*ones(size(bStts))];
    root.burst_len = [root.burst_len; bLens'];
end

[root.burst_tsb,sortinds] = sort(root.burst_tsb);
root.burst_cl = root.burst_cl(sortinds);
root.burst_len = root.burst_len(sortinds);

end