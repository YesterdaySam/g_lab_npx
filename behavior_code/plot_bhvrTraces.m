function [fhandle] = plot_bhvrTraces(sess, laps)
% Plots overlaid, normalized position, velocity, licks, rwd delivery for
% all trials specified in laps with a scale bar
% Position normalized to sess.maxPos, velocity normalized to 99th percentile

ct = 1;
fhandle = figure; hold on; axis off
fhandle.Renderer = 'Painters';

for i = laps
    lapInds = sess.lapstt(i):50:sess.lapend(i);     % Downsample 50x
    xvals = sess.ts(lapInds) - sess.ts(sess.lapstt(i));
    if ct == 1
        % Scale bar is 0.5 seconds by normalization of 1
        plot([sess.ts(sess.lapend(i))+0.5 sess.ts(sess.lapend(i))+1.5] - sess.ts(sess.lapstt(i)), [ct ct], 'k')
        plot([sess.ts(sess.lapend(i))+0.5 sess.ts(sess.lapend(i))+0.5] - sess.ts(sess.lapstt(i)),[ct, ct + 1],'k')
    end
    plot(xvals, sess.pos(lapInds)./sess.maxPos + ct,'k');
    plot(xvals, sess.velshft(lapInds)./45 + ct,'r');    %
    plot(xvals, sess.lck(lapInds)./max(sess.lck)./2 + ct,'b');
    plot(sess.ts(sess.rwdind(i-1)) - sess.ts(sess.lapstt(i)), 1+ct, 'kv')
    ct = ct - 1.25;
end

end