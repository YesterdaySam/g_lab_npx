function [binSpike] = binarizeSpikeTrain(root,unit,sess)
%%% Create a binarized spike train of length sess.ts - 1
% Input
%   root = root object
%   unit = ID of cluster, 0-indexed
%   sess = sess object
%
% Output
%   binSpike = binarized spike train, length sess.ts - 1
%
% Created 10/1/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

spikeTS = sess.ts(root.tsb(root.cl == unit));

binSpike = histcounts(spikeTS,sess.ts);

multicounts = binSpike > 1;
binSpike(multicounts) = 1;  % Ignore multiple spikes in the same time bin

end