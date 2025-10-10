function [spkMap] = binarizeUnits(root,units,sess)
%%% Generate a binarized map of spike times for all units
% Input
%   root = root object
%   units = IDs of units to binarize, 0-indexed
%   sess = sess object
%
% Output
%   spkMap = binarized spike matrix, length sess.ts - 1
%
% Created 10/1/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

nUnits = length(units);
spkMap = zeros(nUnits,length(sess.ts)-1);

for i = 1:nUnits
    cc = units(i);

    spkMap(i,:) = binarizeSpikeTrain(root,cc,sess);
end

spkMap = logical(spkMap);   % To save memory

end