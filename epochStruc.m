function [sess2,root2] = epochStruc(sess, root, indrange)
%% Truncate a sess and root struc based on the input session index range
% Inputs
%   sess = session struct
%   root = root struct
%   indrange = [start stop] indices of the parent sess variable
% Output
%   sess2 = truncated session variable spanning the time specified by indrange
%   root2 = truncated root variable with data within the time of indrange
%
% Created 5/31/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

%% Session
sess2 = sess;
rstInd = indrange(1) - 1;
indrange = indrange(1):indrange(2);

% Analog data
sess2.aidat     = sess.aidat(indrange,:);
sess2.ind       = sess.ind(indrange) - rstInd;
sess2.ts        = sess.ts(indrange);
sess2.vel       = sess.vel(indrange);
sess2.velshft   = sess.velshft(indrange);
sess2.pos       = sess.pos(indrange);
sess2.lck       = sess.lck(indrange);
sess2.lckind    = sess.lckind(sess.lckind > indrange(1) & sess.lckind < indrange(end)) - rstInd;
sess2.opto      = sess.opto(indrange);
sess2.optoind   = sess.optoind(sess.optoind > indrange(1) & sess.optoind < indrange(end))  - rstInd;
sess2.slx       = sess.slx(indrange);
sess2.runInds   = sess.runInds(indrange);
sess2.lapInclude= sess.lapInclude(indrange);

% Digital data
sess2.didat     = sess.didat(indrange);
sess2.rwd       = sess.rwd(indrange);
sess2.rwdind    = sess.rwdind(sess.rwdind >= indrange(1) & sess.rwdind <= indrange(end)) - rstInd;
sess2.rst       = sess.rst(indrange);

% Laps
sess2.lapstt    = sess.lapstt(sess.lapstt >= indrange(1) & sess.lapstt <= indrange(end)) - rstInd;
sess2.lapend    = sess.lapend(sess.lapstt >= indrange(1) & sess.lapstt <= indrange(end)) - rstInd;
sess2.nlaps     = size(sess2.lapend,1);
sess2           = getErrorTrials(sess2);

%% Root
root2 = root;
indroot = [find(root.tsb >= indrange(1),1,'first') find(root.tsb <= indrange(end),1,'last')];
indtimes = sess.ts(root.tsb(indroot));
indlfp  = [find(root.lfp_tsb >= indrange(1),1) find(root.lfp_tsb <= indrange(end),1,'last')];

root2.ts        = root.ts(indroot(1):indroot(2)) - indtimes(1);
root2.cl        = root.cl(indroot(1):indroot(2));
root2.syncpulse = sess2.slx;    % Crude replacement
root2.tspulse   = sess2.ts;     % Crude replacement
root2.lfp       = root.lfp(:,indlfp(1):indlfp(2));
root2.tsb       = root.tsb(indroot(1):indroot(2))  - rstInd;
root2.lfp_tsb   = root.lfp_tsb(indlfp(1):indlfp(2)) - rstInd;
root2.thEnv     = root.thEnv(:,indlfp(1):indlfp(2));
for i = 1:length(root2.ripStruc)
    tmpripples = root.ripStruc(i).ripples;
    root2.ripStruc(i).ripples = tmpripples(tmpripples(:,1) >= root.lfp_tsb(indlfp(1)) & tmpripples(:,3) <= root.lfp_tsb(indlfp(2)),:);
    root2.ripStruc(i).ripples(:,1:3) = root2.ripStruc(i).ripples(:,1:3) - rstInd;
end

end
