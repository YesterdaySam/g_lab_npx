function [lapScore] = get_trialSuccess(sess,useLaps)
%% Calculate 
% Inputs
%   sess    = struct from importBhvr.m
%
% Outputs
%   lapScore = [1xN] of all trials rewarded (1) or not (0)
%
% Created 4/1/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    useLaps = []
end

if isempty(useLaps)
    useLaps = sess.valTrials;
end
nlap = length(useLaps);

for i = 1:nlap
    lapScore(i) = ~isempty(find(sess.rwdTrials == useLaps(i), 1));
end

end