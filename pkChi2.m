function [p,stats,pkBin] = pkChi2(tmpPks,binedges)

binsz = binedges(2)-binedges(1);
binpos = binedges(1:end-1)+binsz/2;
nBins = length(binpos);

for i = 1:size(tmpPks,1)
    pkBin(i) = binpos(logical(tmpPks(i,:)));
end

expec = length(pkBin)/nBins*ones(nBins,1);  % Uniform random counts based on M units and N bins
[~,p,stats] = chi2gof(pkBin,'Expected',expec,'Edges',binedges);

end