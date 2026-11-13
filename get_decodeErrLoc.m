function [errLoc,pErr] = get_decodeErrLoc(dErr,binedges)

arguments
    dErr
    binedges = 0:0.05:1.85
end

errDist = histcounts(abs(dErr),binedges);
binpos = binedges(2:end) - 0.5*diff(binedges(1:2));

[~,pkLoc] = max(errDist);

errLoc = binpos(pkLoc);
pErr = errDist ./ sum(errDist);

end