function [sess2] = epochsess(sess, indrange)
%%
% Inputs
%   indrange = [start stop] indices of the parent sess variable
% Output
%   sess2 = truncated session variable spanning the time specified by indrange
%--------------------------------------------------------------------------

sess2 = sess;
indrange = indrange(1):indrange(2);

% Analog data
sess2.aidat     = sess.aidat(indrange,:);
sess2.ind       = sess.ind(indrange);
sess2.ts        = sess.ts(indrange);
sess2.vel       = sess.vel(indrange);
sess2.pos       = sess.pos(indrange);
sess2.lck       = sess.lck(indrange);
sess2.lckind    = sess.lckind(sess.lckind > indrange(1) & sess.lckind < indrange(end));
sess2.slx       = sess.slx(indrange);

% Digital data
sess2.didat      = sess.didat(indrange);
sess2.rwd        = sess.rwd(indrange);
sess2.rwdind     = sess.rwdind(find(sess.rwdind > indrange(1) & sess.rwdind < indrange(end)));
sess2.rst        = sess.rst(indrange);
sess2.lapstt     = sess.lapstt(find(sess.lapstt > indrange(1) & sess.lapstt < indrange(end)));
sess2.lapend     = sess2.lapstt - 1;
sess2.nlaps      = size(sess2.lapend,1);

end
