function [sess] = get_RunInds(sess,vthresh,tthresh)
%% Identifies running periods and creates sess.runInds subfield
%
% Inputs:
% sess = session struct from importBhvr
% vthresh = threshold of behavioral velocity to throw out, default 0.04 m/s
% tthresh = time threshold of long standing periods to throw out
%
% Outputs:
% sess = modified session struct with new subfield sess.runInds added
%
% Created 12/05/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    vthresh = 0.04;     %m/s
    tthresh = 2;        %s
end

restInds = find(sess.velshft < vthresh);

rests(:,1) = [1; restInds+1];
rests(:,2) = [restInds-1; length(sess.velshft)];

% restDiff = diff(sess.ts(restInds));  %Seconds between stationary periods
% restOns  = find(restDiff > tthresh); 

rests(rests(:,2)-rests(:,1) < 0, :) = [];   %Remove contiguous crossings

rests = reshape(rests', 1, numel(rests)); % reshape into vector like: start, stop, start2, stop2, start3...

tmpdiff = diff(rests);
tmpdiff(1:2:end) = 1+5; %make all inter-epoch times greater than tthresh

mergeinds = find(tmpdiff<1);

mergeinds = cat(2, mergeinds, mergeinds+1);

rests(mergeinds) = [];    % remove all overlapping epochs

rests = reshape(rests, 2, numel(rests)/2)';

rests(rests(:,2) - rests(:,1) < tthresh,:) = [];  % Remove rests smaller than tthresh

tf = zeros(length(sess.velshft),1); % build tf logical vector

for i = 1:size(rests,1)
    
    tf(rests(i,1):rests(i,2)) = 1;
    
end

tf = logical(tf);

sess.runInds = tf;
end