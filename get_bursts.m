function [burstStts, burstEnds, burstLens, burstuISI] = get_bursts(root,sess,unit,tthresh,bthresh)
% Finds bursts of spikes under thresh latency, return start, stop, length
%
% Inputs
%   root = root struct
%   sess = session struct
%   unit = single unit ID
%   plotflag = 0 or 1, whether to plot
%   tthresh = ms, default 8ms between spikes (Bohm et al., 2015)
%   bthresh = # spikes / burst, default 3
%
% Outputs:
%   burstStts = Indices of burst onsets, relative to sess.ind
%   burstEnds = Indices of burst offsets, relative to sess.ind
%   burstLens = number of spikes in burst
%   burstISIs = average ISI of each burst (ms)
%
% Created 12/1/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    sess
    unit
    tthresh = 0.008  % ms
    bthresh = 3
end

spkinds = root.tsb(root.cl == unit);
spkts = sess.ts(spkinds);

spkDiff1 = diff(spkts);
lowLat = double(spkDiff1 < tthresh);

try
    [~,bOnset] = findpeaks(lowLat);
    [~,bOffset] = findpeaks(double(~logical(lowLat))); % Get peaks on inverse array
catch
    disp(['Not enough spikes for burst detection in unit ' num2str(unit)])
    burstStts = NaN;
    burstEnds = NaN;
    burstLens = NaN;
    burstuISI = NaN;
    return
end

if lowLat(1) == 0   % If not starting on a burst
    bLength = (bOffset - bOnset(1:length(bOffset))) + 1;
else    % If starting on a burst
    bOnset = [1 bOnset];
    bLength = (bOffset - bOnset(1:length(bOffset))) + 1; % Add 1 for first spike
end

bursts = bLength >= bthresh; % Trains of 3+ spikes
burstStts = spkinds(bOnset(logical([bursts, 0])));
burstEnds = spkinds(bOffset(logical([bursts, 0])));
burstLens = bLength(bursts);
burstuISI = mean(1000*((sess.ts(burstEnds) - sess.ts(burstStts)) ./ burstLens));

end