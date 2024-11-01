function [root2] = shiftTrain(root,sess)


arguments
    root
    sess
end

dur = sess.ts(end) - sess.ts(1);
nts = size(sess.ts,2);

minShift = round(0.5 * sess.samprate);     %Minimum duration to shift spike trains (e.g. 500ms)
rShift = minShift + randi(nts - 2*minShift,1);    % Get # of tsb to shift all spike trains, but don't over-shift

root2 = root;
root2.tsb = root2.tsb + rShift;         % Shift spikes
wrapSpks = root2.tsb > sess.ind(end);    % Find t-shifted spikes greater than session endpoint

root2.tsb(wrapSpks) = root2.tsb(wrapSpks) - sess.ind(end);

end