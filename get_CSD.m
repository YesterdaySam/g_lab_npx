function [] = get_CSD(root,band,refch)
% Calculates second spatial derivative of LFP in band across depth
% Based on Ylinen et al., 1995 and Bragin et al., 1995
% Not 100% certain this is working - need to check how spatial derivative
% is taken
% Inputs
%   root = root object, must have multiple channels worth of LFP (i.e. use addRootLFP)
%   band = [low high] vector of frequencies e.g. [6 10]
% 
% Outputs
%
% Created 1/20/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

rng(1);

filtlfp = bandpass(root.lfp', band, root.fs_lfp);
sigEnv = abs(hilbert(filtlfp(:,refch)));

twin = round(root.fs_lfp/4);

    % find peaks of band and sample 100x 'trials'
[pks,locs] = findpeaks(sigEnv);
hiPks = pks > prctile(pks,90);  % Get peaks over 90th percentile
hiLocs = locs(hiPks); edgeLocs = hiLocs < twin & hiLocs > root.lfp(end) - twin; hiLocs(edgeLocs) = [];
locDiff = diff(hiLocs);
distLocs = locDiff > root.fs_lfp;
randsamps = randperm(length(distLocs),100); % No replacement
randLocs = sort(hiLocs(randsamps));

lfMap = [];
for i = 1:length(randLocs)
    lfMap(:,:,i) = filtlfp(randLocs(i)-twin:randLocs(i)+twin,:);
end
lfMean = mean(lfMap,3); % Average over the 'trials'

depthSmooth = movmean(lfMean,3,2);
dy2 = -diff(depthSmooth,2,2);    % Take second derivative

figure;
imagesc(depthSmooth'); set(gca,'YDir','normal'); colormap(turbo)
title('Smoothed, Avg LFP')
figure;
imagesc(dy2'); set(gca,'YDir','normal'); colormap(turbo)
title('CSD')
end