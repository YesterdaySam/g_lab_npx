function [lapScore] = get_trialSuccessVshuffle(sess,rzloc,useLaps,nshufs)
%% Calculate 
% Inputs
%   sess    = struct from importBhvr.m
%   rzloc   = in m, position of reward zone start
%   useLaps = [] or array of laps to use, default is all sess.valTrials
%   nshufs  = 250; number of shuffles to perform
%
% Outputs
%   lapScore = Real - chance lap reward perforamnce
%
% Created 10/14/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    rzloc
    useLaps = []
    nshufs = 250
end

rng(1) % Fix rng seed

realSucc = get_trialSuccess(sess,useLaps);
uRealSucc = mean(realSucc,'omitnan');

if isempty(useLaps)
    useLaps = sess.valTrials;
end
nlap = length(useLaps);

zoneSize = 0.3;
zonePass = 0;
rzEnd = rzloc + zoneSize;

for i = 1:nshufs
    for j = 1:nlap
        shufZone = [];
        while zonePass == 0
            shufStt = (sess.maxPos) * rand(1); % Pick a random zone start
            % Verify if the shuffle zone start is outside of a 90cm window from the
            % real zone end, representing 30 * 3 cm of real RZ, anticipatory, and 30cm shuffle zone end
            if (rzEnd - shufStt > 0.9 || rzEnd - shufStt < 0) && ~(rzEnd - shufStt < -0.95)
                zonePass = 1;
                % tmpStt(j) = shufStt; % Storage for testing zone starts
            end
        end
        zonePass = 0;

        shufEnd = shufStt + zoneSize;
        if shufEnd > sess.maxPos
            wrapZone = [0 zoneSize - (sess.maxPos - shufStt)];
            shufZone(1,:) = wrapZone;
            shufZone(2,:) = [shufStt sess.maxPos];
            % tmpEnd(j) = wrapZone(2);  % Storage for testing zone starts
        else
            shufZone = [shufStt shufStt + zoneSize];
            % tmpEnd(j) = shufStt + zoneSize;  % Storage for testing zone starts
        end
        
        lapLcks = sess.lckind(sess.lckind > sess.lapstt(useLaps(j)) & sess.lckind < sess.lapend(useLaps(j)));
        for k = 1:size(shufZone,1)
            lcksIn(k) = histcounts(sess.pos(lapLcks), shufZone(k,:));
        end

        shufSuccess(j) = sum(lcksIn) > 0;
    end
    uShufSucc(i) = mean(shufSuccess,'omitnan');
end

% Sanity check 1 shuffle worth of start/stop times
% figure; histogram(tmpStt,0:0.05:1.85)
% hold on; histogram(tmpEnd,0:0.05:1.85)

lapScore = uRealSucc - mean(uShufSucc);

end