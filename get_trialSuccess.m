function [lapScore] = get_trialSuccess(sess)
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
end

for i = 1:length(sess.valTrials)
    lapScore(i) = ~isempty(find(sess.rwdTrials == sess.valTrials(i), 1));
end

end