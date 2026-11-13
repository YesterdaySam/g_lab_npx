function [alnSpSpwrRate,alnSpSpwr] = get_spatialSPWRs(root,sess,datstruc,rzloc)
%% Bin SPWRs by spatial bins in datstruc, align to rzloc
% Inputs:
%   root = root struc
%   sess = sess struc
%   datstruc = data with shuffles
%   rzloc = distance (m) of reward zone, default is not to shift
%
% Outputs:
%   alnSpSpwr = probability normalized by all spatial bins, of ripples
%   alnSpSpwrRate = SPWR rate by spatial position in datstruc.binedges
%
% Updated 10/27/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    sess
    datstruc
    rzloc = round(length(datstruc.binpos)/2)   % Default don't shift
end

mBin = round(length(datstruc.binpos)/2);

rzbin = find(datstruc.binpos > rzloc,1);
spSpwrs = histcounts(sess.pos(root.ripStruc(root.ripRef).ripples(:,2)),datstruc.binedges); % Get rip count
spStand = histcounts(sess.pos(~sess.runInds),datstruc.binedges)/sess.samprate; % Get standing time (secs)
alnSpSpwrRate = circshift(spSpwrs./spStand, mBin - rzbin);
alnSpSpwrRate(isnan(alnSpSpwrRate)) = 0;
normSpwrs = spSpwrs ./ sum(spSpwrs);
alnSpSpwr = circshift(normSpwrs, mBin - rzbin);

end