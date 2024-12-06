function [root] = shiftTrain(root,sess)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% sess = session struct from importBhvr
%
% Outputs:
% root = modified root object with circularly shifted root.tsb field
% Created 9/19/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------


arguments
    root
    sess
end

% dur = sess.ind(end) - sess.ind(1);
dur = length(root.tspulse);
nts = size(sess.ts,2);

minShift = round(0.5 * sess.samprate);     %Minimum duration to shift spike trains (e.g. 500ms)
rShift = minShift + randi(nts - 2*minShift,1);    % Get # of tsb to shift all spike trains, but don't over-shift

% Shift spikes
root.tsb = mod(root.tsb + rShift, dur-1)+1;         % Use Modulo for increased speed, and offset by 1 to avoid indexing '0's

% root.tsb = root.tsb + rShift;         % Shift spikes
% wrapSpks = root.tsb > sess.ind(end);    % Find t-shifted spikes greater than session endpoint
% 
% root.tsb(wrapSpks) = root.tsb(wrapSpks) - sess.ind(end);

end