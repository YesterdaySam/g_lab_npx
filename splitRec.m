function [sessPre, sessPst, rootPre, rootPst] = splitRec(sess,root,splitLap,useLap)
%% 
% Inputs:
%   root = root object
%   sess = sess object
%   splitLap = optional; can specify lap number
%
% Outputs:
%   sessPre = session object epoched before splitLap
%   sessPst = session object epoched after  splitLap
%   rootPre = session object epoched before splitLap
%   rootPst = session object epoched after  splitLap
%
% Created 6/4/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    root = 0
    splitLap = 0
    useLap = []
end

if splitLap == 0
    doRMNonOps = 0;
    splitLap = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
    splitLap = sess.rwdTrials(splitLap);
else
    doRMNonOps = 1;
    nonops     = get_nonOpLaps(sess);
    sess       = get_errorTrials(sess,nonops);
    sess       = get_lapInclude(sess);
end

if isempty(useLap)
    lapThreshF = sess.ind(1);
    lapThreshN = sess.ind(end);
else
    valF = sess.valTrials(sess.valTrials <= splitLap);
    valN = sess.valTrials(sess.valTrials > splitLap);
    lapThreshF = sess.lapstt(valF(end-useLap+1));
    lapThreshN = sess.lapend(valN(useLap));
end

frstHalfInds       = [lapThreshF sess.lapend(splitLap)];
lastHalfInds       = [sess.lapstt(splitLap+1) lapThreshN];
[sessPre, rootPre] = epochStruc(sess,root,frstHalfInds,doRMNonOps);
[sessPst, rootPst] = epochStruc(sess,root,lastHalfInds,doRMNonOps);

figure; hold on; 
set(gcf,'Position',[1302 823 1450 420])
plot(sess.pos)
plot(sess.ind(frstHalfInds(1):frstHalfInds(2)),sessPre.pos)
plot(sess.ind(lastHalfInds(1):lastHalfInds(2)),sessPst.pos)
fsave(gcf,'newsplit',1)
close all
end