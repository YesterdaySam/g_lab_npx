function [] = get_bursts(root,sess,cc,thresh)

spkts = sess.ts(root.tsb(root.cl(cc)));

spkDiff1 = diff(spkts);
lowLat = spkDiff1 < thresh;

spkDiff2 = diff
end