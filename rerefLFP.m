function [root] = rerefLFP(root,lfp,depths,shanks)
% shanks = Use 1-4 for shanks

for i = 1:length(shanks)
    chan = find(lfp.lfpinfo.lfpDepth == depths(i) & lfp.lfpinfo.lfpShank == shanks(i)-1);
    root.lfp(shanks(i),:) = lfp.lfp(chan,:);
end

end