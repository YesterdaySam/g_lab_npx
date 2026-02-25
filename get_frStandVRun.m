function [binfr] = get_frStandVRun(root,unit,sess,vthresh)
%% Calculate mean firing rate of [standing, running] periods
%
% Inputs:
% root = root object. Must have uPSDMax field
% unit = cluster ID
% sess = session struct from importBhvr
% vbnsz = size of velocity bins, default 2cm/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% binfr = occupancy-normalized firing rate of standing and running periods
%
% Created 3/15/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    vthresh = 4       %cm/s
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);

binedges = [-vthresh, vthresh, max(sess.velshft)];

spkvel = sess.velshft(root.tsb(root.cl == unit));

vspk = histcounts(spkvel,binedges);
vocc = histcounts(sess.velshft,binedges)/sess.samprate;
binfr = vspk ./ vocc;

end