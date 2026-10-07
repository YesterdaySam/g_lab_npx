function [nonOpLaps] = get_nonOpLaps(sess)

minLick = [];
for j = 1:length(sess.rwdind)
    minLick(j) = min(abs(sess.rwdind(j) - sess.lckind));
end

nonOpLaps = sess.rwdTrials(minLick > 10);

end