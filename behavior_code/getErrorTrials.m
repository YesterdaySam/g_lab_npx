function [sess] = getErrorTrials(sess,manualtrials)
% Get error, valid, and rewarded trials and assign them to sess variable
%
% Inputs
%   sess = from importBhvr.m
%   manualTrials = array of indices of 
%
% Outputs
%   sess = modified with sess.errTrials (long or short)
%                        sess.valTrials (right length only)
%                        sess.rwdTrials (right length & rewarded)
%
% Created 8/6/24 LKW; Grienberger Lab; Brandeis University
% Modified 10/22/25 LKW
%--------------------------------------------------------------------------

arguments
    sess                %session struct
    manualtrials = []   %Nx1 array of user-defined bad trials
end

if isfield(sess,'maxPos')
    maxPos = sess.maxPos;
else
    maxPos = 1.86;
end

isgood = ones(size(sess.lapstt));

% Find long or short laps and exclude
for i = 1:length(sess.lapstt)
    lapLen(i) = max(sess.pos(sess.lapstt(i):sess.lapend(i))) - min(sess.pos(sess.lapstt(i):sess.lapend(i)));

    if lapLen(i) > maxPos + 0.2    % in meters
        isgood(i) = 0;
    elseif lapLen(i) < maxPos - 0.2
        isgood(i) = 0;
    end
end

% Exclude user-defined laps
isgood(manualtrials) = 0;

% Find rewarded laps, starting with ts(1) and ending on ts(end)
isrwd = logical(histcounts(sess.rwdind,[sess.lapstt(1); sess.lapend]));

% Assign to sess
isgood = isgood > 0;
trlist = 1:sess.nlaps;
sess.errTrials = trlist(logical(1-isgood));
sess.valTrials = trlist(logical(isgood));
sess.rwdTrials = trlist(isgood & isrwd');
end