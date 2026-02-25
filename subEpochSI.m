function [datStruc] = subEpochSI(sess,root,datStruc,nPer)
% Breaks sess into groups of nPer laps and returns Spatial Info for each
% Inputs
%   sess = sess object
%   root = root object
%   datStruc = data structure to save resulting SI
%   nPer = integer number of laps over which to calculate SI
% Outputs
%   datStruc = modified data structure with SI

arguments
    sess
    root
    datStruc
    nPer = 10   % Groups of 10 laps
end

dbnsz = 0.05;

nLaps = length(sess.valTrials);
epochSubInds = [sess.lapstt(sess.valTrials(1:nPer:nLaps-mod(nLaps,nPer))),...
    sess.lapend(sess.valTrials(1:nPer:nLaps-mod(nLaps,nPer))+(nPer-1))];

for i = 1:size(epochSubInds,1)
    [sesstmp, roottmp] = epochStruc(sess,root,epochSubInds(i,:));
    for j = 1:length(root.good)
        cc = root.good(j);
        [datStruc.subEpochSI(j,i)] = get_SI(roottmp,cc,sesstmp,dbnsz);
    end
end

end