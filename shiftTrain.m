function [root,sess,rShift] = shiftTrain(root,sess)
%% Circularly shifts the tsb field of the root object 
%
% Inputs:
% root = root object.
% sess = session struct from importBhvr
%
% Outputs:
% root = modified root object with circularly shifted root.tsb field
% sess = modified sess object -- not validated!
% rShift = index of pseudo-random shift
%
% Created 9/19/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------


arguments
    root
    sess
end

% dur = sess.ind(end) - sess.ind(1);
dur = length(root.syncpulse);
nts = size(sess.ts,2);

minShift = round(0.5 * sess.samprate);     %Minimum duration to shift spike trains (e.g. 500ms)
rShift = minShift + randi(nts - 2*minShift,1);    % Get # of tsb to shift all spike trains, but don't over-shift

% Shift spikes
root.tsb = mod(root.tsb + rShift, dur-1)+1;         % Use Modulo for increased speed, and offset by 1 to avoid indexing '0's
sess.runInds = [sess.runInds(end-rShift:end); sess.runInds(1:end-rShift-1)];

% root.tsb = root.tsb + rShift;         % Shift spikes
% wrapSpks = root.tsb > sess.ind(end);    % Find t-shifted spikes greater than session endpoint
% 
% root.tsb(wrapSpks) = root.tsb(wrapSpks) - sess.ind(end);


%%
% dur = root.tsb(end);
% % dur = length(root.tspulse);
% nts = size(sess.ts,2);
% 
% minShift = round(0.5 * sess.samprate);     % Minimum duration to shift spike trains (e.g. 500ms)
% rShift = minShift + randi(nts - 2*minShift,1);    % Get # of tsb to shift all spike trains, but don't over-shift
% % rShift = 100;
% % tmptsb = root.tsb;
% 
% % Shift spikes
% root.tsb = mod(root.tsb + rShift, dur + rShift) +1;         % Use Modulo for increased speed, and offset by 1 to avoid indexing '0's
% sess.runInds = [sess.runInds(end-rShift:end); sess.runInds(1:end-rShift-1)];
% 
% % root.tsb = root.tsb + rShift;         % Shift spikes
% % wrapSpks = root.tsb > sess.ind(end);    % Find t-shifted spikes greater than session endpoint
% % 
% % root.tsb(wrapSpks) = root.tsb(wrapSpks) - sess.ind(end);

end