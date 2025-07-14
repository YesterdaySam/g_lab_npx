function [vMap] = get_periRZVel(sess,rzPos,rzLen)
%% Calculate Velocity in peri-RZ location - Not validated!
% Inputs
%   sess    = struct from importBhvr.m
%   rzPos   = [1x2] in cm of [RZ position, comparison position]
%   rzLen   = double of pre-RZ length in cm
%
% Outputs
%   vMap    = peri-RZ avg velocity by trial
%
% Created 7/8/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    rzPos           % in cm, [1x2] of rZ positions
    rzLen = 30      % in cm
end

rzPos = rzPos ./ 100;
rzLen = rzLen ./ 100;

rzPosMat = [rzPos(1) - rzLen, rzPos(1)];
if rzPosMat(1) < 0   % if RZ length puts the start of RZ in previous lap
    rzPosMat(1,:) = [0, rzPos(1)];
    rzPosMat(2,:) = [max(sess.pos(sess.lapInclude)) + (rzPos(1) - rzLen), max(sess.pos(sess.lapInclude))];
end
rzPosMat = [rzPosMat; rzPos(2) - rzLen, rzPos(2)];

[binedges,bnvel] = plot_trialvel(sess,0.01,0);

for i = 1:size(rzPosMat,1)
    [~,loc(i,:)] = min(abs(rzPosMat(i,:) - binedges'));
end
% [~,loc(2,:)] = min(abs(rzPosMat(2,:) - binedges'));
% [~,loc(3,:)] = min(abs(rz2(1,:) - binedges'));

for i = 1:size(loc,1)
    velMap(:,i) = mean(bnvel(:,loc(i,1):loc(i,2)-1),2,'omitnan');
end

if size(velMap,2) == 3
    vMap(:,1) = mean([velMap(:,1), [velMap(1,2); velMap(1:end-1,2)]],2);
    vMap(:,2) = velMap(:,3);
elseif size(velMap,2) == 2
    vMap = velMap;
else
    disp('error, more than 3 dimensions to velMap')
    vMap = [];
end

end