function [root] = get_lfpXdepth(root, bands, trange)
%% Calculates mean Welch's Power Spectral Density (PSD) for certain bands
%
% Inputs:
% root = root object. Must have root.thEnv field (bandpassed, hilbert
%       transformed theta)
% bands = Nx2 double of [freq_lower, freq_higher] where N is number of
%       bandpass frequency ranges over which to estimate power
% trange = indices of root.fs_lfp over which to estimate; e.g. [start_ind, stop_ind]
%       Useful if artifactual noise is preventing a clean PSD estimation
%
% Outputs:
% root = updated root object including the following fields:
%       bands = Nx2 double of [freq_lower, freq_higher] used in PSD estimate
%       uPSD = mean PSD for band in bands for each LFP channel in root.lfp
%       uPSDMax = max PSD for band in bands by probe shank
%       thEnv = envelope of theta band by shank, assuming bands(1,:) includes theta
%
% Created 5/8/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    bands = [6 10; 150 250];    % Spectral bands to estimate; Also can try [0.01 300]
    trange = [];    % In lfp_tsb indices. If user supplied, then only estimate power over that time range
end

if ~isempty(trange)
    tmpLFP = root.lfp(:,trange(1):trange(2));
else
    tmpLFP = root.lfp;
end

nshanks = numel(unique(root.info.shankID));
NFFT = 10000;
[pxx,f] = pwelch(tmpLFP',[],[],NFFT,root.fs_lfp);

uPSD = zeros(size(bands,1),size(pxx,2));

for band = 1:size(bands,1)
    f_tmp = f > bands(band,1) & f < bands(band,2);  % Get frequencies in band
    uPSD(band,:) = mean(10*log10(pxx(f_tmp,:)),1);    %Get mean per electrode within band
end

% Locate electrode with max spectral power for each band, for each shank
root.bands = bands;
for sh = 1:nshanks
    for band = 1:size(bands,1)
        [tmpMax, tmpInd] = max(uPSD(band,root.lfpinfo.lfpShank == sh-1));
        tmpD = root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh-1);
        tmpD = tmpD(tmpInd);
        root.uPSDMax(band,sh) = find(root.lfpinfo.lfpDepth == tmpD & root.lfpinfo.lfpShank == sh-1);
    end
end

root.uPSD = uPSD;

% Interpolate over bad channels with low PSD 
badChans = root.uPSD(2,:) < mean(root.uPSD(2,:)) - 5* std(root.uPSD(2,:)); % With average LFP < 5 StDev of mean
if sum(badChans) > 0
    disp(['Found ' num2str(sum(badChans)) ' Bad LFP channels, spline interpolation commencing'])
    for sh = 1:nshanks   %Interpolate by shank instead of globally
        shChans = root.lfpinfo.lfpShank == sh-1;
        tmpsh_upsd = root.uPSD(2,shChans);
        badChans = tmpsh_upsd < mean(tmpsh_upsd) - 3* std(tmpsh_upsd); % With average LFP < 3 StDev of mean
        if sum(badChans) > 0 
            interpVec = 1:length(tmpsh_upsd);
            interpTrn = interpVec(~badChans);
            uPSDInterp = interp1(interpTrn,tmpsh_upsd(~badChans),interpVec,'spline');
            root.uPSD(2,shChans) = uPSDInterp;
        end
    end
end

% Pre-calculate theta at s.pyr. center, based on 150-250Hz PSD peak across shanks
for i = 1:nshanks
    thetaLFPSh = bandpass(root.lfp(root.uPSDMax(2,i),:),root.bands(1,:), root.fs_lfp);
    root.thEnv(i,:) = hilbert(thetaLFPSh);  % Operates along rows
end

end