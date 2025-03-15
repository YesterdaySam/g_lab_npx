function [nearInds] = compare_ripples_times(t1s,t2s,fs,tthresh)
%% Collects nearest time point of t2 to each t1 time under tthresh
%
% Inputs:
% t1s = indices of event times for first group
% t2s = indices of event times for second group to be compared to first
% fs  = sampling rate of event times
% tthresh = double in seconds; cutoff for t2s too far from t1s
%
% Outputs:
% nearInds = vector of length t1s where each element is the nearest index
% from t2s, or NaN if no element closer than tthresh*fs
%
% Created 3/5/2025 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    t1s {double}   % Vector of time indices for reference
    t2s {double}   % Vector of time indices for comparison
    fs      = 1250 % sampling rate
    tthresh = 1    % 1sec
end

nEvents1 = length(t1s);

nearInds = zeros(nEvents1,1);
for i = 1:nEvents1
    [tmpMin, tmpMinInd] = min(abs(t1s(i) - t2s));
    if tmpMin < tthresh*fs  % If under time threshold
        nearInds(i) = t2s(tmpMinInd);
    else
        nearInds(i) = NaN;
    end
end

end