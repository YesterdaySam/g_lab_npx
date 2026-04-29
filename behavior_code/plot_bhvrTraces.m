function [fhandle] = plot_bhvrTraces(sess, laps)
% Plots overlaid, normalized position, velocity, licks, rwd delivery for
% all trials specified in laps with a scale bar
% Position normalized to sess.maxPos, velocity normalized to 99th percentile

ctpos = 1;
ctlap = 1;
fhandle = figure; hold on; axis off
fhandle.Renderer = 'Painters';

for i = laps
    if ctlap > 1    % Skip a Y-position if there's a gap between laps
        if laps(ctlap) - laps(ctlap-1) > 1
            ctpos = ctpos - 1.25;
        end
    end

    lapInds = sess.lapstt(i):50:sess.lapend(i);     % Downsample 50x
    xvals = sess.ts(lapInds) - sess.ts(sess.lapstt(i));
    lcks = sess.lckind(sess.lckind > sess.lapstt(i) & sess.lckind < sess.lapend(i));
    if ctpos == 1
        % Scale bar is 0.5 seconds by normalization of 1
        plot([sess.ts(sess.lapend(i))+0.5 sess.ts(sess.lapend(i))+1.5] - sess.ts(sess.lapstt(i)), [ctpos ctpos], 'k')
        plot([sess.ts(sess.lapend(i))+0.5 sess.ts(sess.lapend(i))+0.5] - sess.ts(sess.lapstt(i)),[ctpos, ctpos + 1],'k')
    end
    plot(xvals, sess.pos(lapInds)./sess.maxPos + ctpos,'k');
    plot(xvals, sess.velshft(lapInds)./45 + ctpos,'r');    %
    % plot(xvals, sess.lck(lapInds)./max(sess.lck)./2 + ct,'b');
    plot(sess.ts(lcks) - sess.ts(sess.lapstt(i)), ones(size(lcks)) + ctpos-0.2, 'b|');
    plot(sess.ts(sess.rwdind(i-1)) - sess.ts(sess.lapstt(i)), 1+ctpos, 'kv')
    ctpos = ctpos - 1.25;
    ctlap = ctlap + 1;
end

end