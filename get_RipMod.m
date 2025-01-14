function [ripMod,uRipFR,uCtlFR,ripSpkMap] = get_RipMod(root,unit,tOffset)
%% Returns the Ripple Modulation Index of a unit
%
% Gets ripple_spks / (ripple_spks + control_spks)
% Where control_spks is a matched duration shifted by tOffset from each ripple
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% tOffset = distance of offset in msec, can be - or +
%
% Outputs:
% ripMod = ripple_spks / (ripple_spks + control_spks)
% uRipFR = firing rate spks/sec average within ripples
% uCtlFR = firing rate spks/sec average within control periods
% ripSpkMap = [n x 2] array of spike # for [ripple control] periods
%
% Created 1/14/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

lowSpkThresh = 25;

ripStt = root.ripples(:,1);
ripEnd = root.ripples(:,3);

offInd = tOffset / 1000 * root.fs_lfp;

ctlStt = ripStt + offInd;
ctlEnd = ripEnd + offInd;

% Find and remove ripples with control periods falling before or after session
badInds = ctlStt < root.lfp_tsb(1) | ctlStt > root.lfp_tsb(end);
badInds = logical(badInds + ctlEnd < root.lfp_tsb(1) | ctlEnd > root.lfp_tsb(end));

ctlStt(badInds) = [];
ctlEnd(badInds) = [];
ripStt(badInds) = [];
ripEnd(badInds) = [];

nRips = length(ripEnd);
spkinds = root.tsb(root.cl == unit);
ripDurs = (ripEnd - ripStt) / root.fs_lfp;    % In seconds

ripSpkMap = zeros(nRips,2);
ripFRMap = zeros(nRips,2);

for i = 1:nRips
    ripSpkMap(i,1) = sum(spkinds > ripStt(i) & spkinds < ripEnd(i));
    ripSpkMap(i,2) = sum(spkinds > ctlStt(i) & spkinds < ctlEnd(i));
    ripFRMap(i,:) = ripSpkMap(i,:) ./ ripDurs(i);
end

if sum(ripSpkMap,'all') < lowSpkThresh
    ripMod = NaN;
    uRipFR = NaN;
    uCtlFR = NaN;
    disp(['Fewer than ' num2str(lowSpkThresh) ' spikes found for unit ' num2str(unit) ' during get_RipMod'])
else
    ripMod = sum(ripSpkMap(:,1)) ./ sum(ripSpkMap,'all');
    uRipFR = mean(ripFRMap(:,1));
    uCtlFR = mean(ripFRMap(:,2));
end

end