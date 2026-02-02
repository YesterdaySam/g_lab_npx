function [bIndex] = get_burstIndex(root,sess,unit)
% Finds bursts of spikes under thresh latency, return start, stop, length
%
% Inputs
%   root = root struct
%   sess = session struct
%   unit = single unit ID
%
% Outputs:
%   bIndex = ratio of burst spikes over non-burst (singlet/doublet) spikes
%
% Created 12/2/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root {struct}
    sess {struct}
    unit {double}
end

[~,~,bLengths] = get_bursts(root,sess,unit);

nBurstSpikes = sum(bLengths);
nUnitSpikes  = sum(root.cl == unit);
nNonBurstSpikes = nUnitSpikes - nBurstSpikes;

bIndex = nBurstSpikes / nNonBurstSpikes;

end