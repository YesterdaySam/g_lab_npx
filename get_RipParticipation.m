function [ripParticipation] = get_RipParticipation(root,unit,sess,wlen)
%% Returns the Ripple Participation Probability of a unit
%
% Counts N ripples where unit fired within wlen, divide by total # ripples
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID OR spike train indices (i.e. when jitter shuffling)
% tOffset = distance of offset in msec, can be - or +
%
% Outputs:
% ripMod = ripple_spks / (ripple_spks + control_spks)
% uRipFR = firing rate spks/sec average within ripples
% uCtlFR = firing rate spks/sec average within control periods
% ripSpkMap = [n x 2] array of spike # for [ripple control] periods
%
% Created 5/14/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit            %Cluster ID OR spike train indices (i.e. when jitter shuffling)
    sess
    wlen     = 25  %msec
end

if isscalar(unit)
    spkinds = root.tsb(root.cl == unit);
    tmpind = find(root.info.cluster_id == unit);
else
    spkinds = unit;
end
spkinds = sess.ts(spkinds);
tmpsz = size(spkinds);
if tmpsz(2) > 1
    spkinds = spkinds';
end

nRips = size(root.ripples,1);
wdw = wlen/1000;
ripSpiked = zeros(nRips, 1);

for i = 1:nRips
    tmpspks      = spkinds(spkinds > sess.ts(root.ripples(i,2)) - wdw & spkinds < sess.ts(root.ripples(i,2)) + wdw) - sess.ts(root.ripples(i,2));
    ripSpiked(i) = ~isempty(tmpspks);
end

ripParticipation = sum(ripSpiked) / nRips;

end