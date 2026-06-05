function [sessPre, sessPst, rootPre, rootPst] = splitRec(sess,root,splitLap)
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
end

if splitLap == 0
    splitLap = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
    splitLap = sess.rwdTrials(splitLap);
end

frstHalfInds       = [sess.ind(1) sess.lapend(splitLap)];
lastHalfInds       = [sess.lapstt(splitLap+1) sess.ind(end)];
[sessPre, rootPre] = epochStruc(sess,root,frstHalfInds);
[sessPst, rootPst] = epochStruc(sess,root,lastHalfInds);

end