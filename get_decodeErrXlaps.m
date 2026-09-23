function [uLapAbsErr,uLapRawErr, errLapLoc] = get_decodeErrXlaps(dcI,sess,binedges)
%% Uses dcI decoder struct to calculate raw and abs error over laps in sess
%
% Inputs:
%   dcI = decoder info struct from decodePosBayes.m
%   sess = sess struct
%   binedges = positional bin vector e.g. 0:0.05:0.95 
%
% Outputs:
%   uLapAbsErr = avg error by lap, absolute distance between real and decode (m)
%   uLapRawErr = avg error by lap, raw distance between real and decode (m)
%   errLapLoc  = bin from binedges with highest density of error distances (m)
%
% Created 9/15/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    dcI
    sess
    binedges = 0:0.05:0.95
end

nLaps = length(sess.lapstt);
binpos = binedges(2:end) - 0.5*diff(binedges(1:2));

uLapAbsErr = nan(nLaps,1);
uLapRawErr = nan(nLaps,1);
errLapLoc     = nan(nLaps,1);

for i = 1:nLaps
    indstt = find(dcI.newT >= sess.ts(sess.lapstt(i)),1);
    indend = find(dcI.newT >= sess.ts(sess.lapend(i)),1) - 1;

    uLapAbsErr(i) = mean(abs(dcI.dErr(indstt:indend)),'omitnan');
    uLapRawErr(i) = mean(dcI.dErr(indstt:indend),'omitmissing');

    errDist = histcounts(abs(dcI.dErr(indstt:indend)),binedges);
    [~,errLapLoc(i)] = max(errDist);

end
errLapLoc = binpos(errLapLoc); % Convert to m from bin inds

end