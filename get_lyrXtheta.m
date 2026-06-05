function [root] = get_lyrXtheta(root,thStats,optoPkT)
%% Assign layer of EC using theta phase preference alone
% Inputs:
%   datT = a table organizing sessions by mouse and recording day
%   fname = save name of the file containing data variables
%   sdir = location to save combined data variables
% Outputs:
%   None (variables saved in-function)
%
% Updated 5/18/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

tmpD = root.info.depth(root.goodind);
[thAng, ~, thSigs] = get_thAng(thStats);
thSig = thSigs < 0.05;

ec2inds = thAng > -90 & thAng < 90;
ec3inds = ~ec2inds; %-180:-90 and 90:180;

dBins = 0:50:max(root.info.depth);
lyrct = histcounts(tmpD,dBins); % All good by depth bin
lyrsg = histcounts(tmpD(thSig),dBins); % Significantly theta by depth
lyrtr2 = histcounts(tmpD(ec2inds),dBins); % Theta tuning by depth
ec5stt = dBins(find(smooth(lyrsg./lyrct) <= 0.1,1)); % Estimate L5 start using cutoff of ratio of significant / all
ec3stt = dBins(find(smooth(lyrtr2./lyrct) <= 0.4,1)); % Estimate L3 start using cutoff of ratio of L2 theta / all

ec5inds = tmpD' > ec5stt & ~thSig;
ec2inds = ec2inds & tmpD' < ec3stt + 300;
try
    ec3inds = (~ec2inds & ~ec5inds) | optoPkT < 0.015;
catch
    ec3inds = ~ec2inds & ~ec5inds;
end

tmpLyr = zeros(size(root.good));
tmpLyr(ec2inds) = 2;
tmpLyr(ec3inds) = 3;
tmpLyr(ec5inds) = 5;
root.info.lyrID(root.goodind) = tmpLyr;

end