function [sess] = getErrorTrials(sess,manualtrials)

arguments
    sess                %session struct
    manualtrials = []   %Nx1 array of user-defined bad trials
end

% Find unrewarded laps, starting with ts(1) and ending on ts(end)
isgood = histcounts(sess.rwdind,[sess.lapstt(1); sess.lapend]);

% Find long laps from failed reset and exclude
for i = 1:length(sess.lapstt)
    trlen = max(sess.pos(sess.lapstt(i):sess.lapend(i)));
    if trlen > 2    % in meters
        isgood(i) = 0;
    end
end

isgood = isgood > 0;
trlist = 1:sess.nlaps;
sess.errTrials = trlist(logical(1-isgood));
sess.valTrials = trlist(logical(isgood));

end