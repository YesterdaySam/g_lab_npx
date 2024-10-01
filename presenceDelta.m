function [presPass, delt] = presenceDelta(root,sess,cc,nBins)

arguments
    root
    sess
    cc
    nBins = 10
    % cutoff = 3
end

spkts = sess.ts(root.tsb(root.cl == cc));

binedges = linspace(sess.ts(1),sess.ts(end),nBins);

bnspks = histcounts(spkts,binedges);

ccMean = mean(bnspks);
ccTotal = sum(bnspks);
ccStdv = std(bnspks);
ctrBin = floor(nBins/2);

delt.FirstLast = diff([bnspks(1) bnspks(end)])/ccTotal;
delt.FirstMid = diff([bnspks(1) bnspks(ctrBin)])/ccTotal;
delt.LastMid = diff([bnspks(ctrBin) bnspks(end)])/ccTotal;

% binZ = (bnspks - ccMean) ./ ccStdv;
% binZz = abs(binZ);

cutoff = 1/nBins;

if abs(delt.FirstLast) > cutoff || abs(delt.FirstMid) > cutoff || abs(delt.LastMid) > cutoff
    presPass = 0;
else
    presPass = 1;
end

end
