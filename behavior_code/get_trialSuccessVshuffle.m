function [lapScore] = get_trialSuccessVshuffle(sess,rzloc,nshufs)
%% Calculate 
% Inputs
%   sess    = struct from importBhvr.m
%
% Outputs
%   lapScore = [1xN] of all trials rewarded (1) or not (0)
%
% Created 10/14/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    rzloc
    nshufs = 250
end

realSucc = get_trialSuccess(sess);
uRealSucc = mean(realSucc,'omitnan');

nind = size(sess.ind,2);
dur = length(sess.ind);
nlap = length(sess.valTrials);

for i = 1:nshufs
    for j = 1:nlap
        shufZone = (sess.maxPos-0.3) * rand(1);
        shufZone = [shufZone shufZone+0.3];

        lapLcks = sess.lckind(sess.lckind > sess.lapstt(sess.valTrials(j)) & sess.lckind < sess.lapend(sess.valTrials(j)));
        lcksIn = histcounts(sess.pos(lapLcks), shufZone);

        shufSuccess(j) = lcksIn > 0;
    end
    uShufSucc(i) = mean(shufSuccess,'omitnan');
end

lapScore = uRealSucc - mean(uShufSucc);

end