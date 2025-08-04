function [lckmap, lickDI] = get_lickDiscrim(sess,rzPos,rzLen)
%% Calculate Discrimiation Index of licks in rzPos(1,:) vs rzPos(2,:)
% Inputs
%   sess    = struct from importBhvr.m
%   rzPos   = [1x2] in cm of [RZ position, comparison position]
%   rzLen   = double of peri-RZ length in cm
%
% Outputs
%   lckmap  = lick counts by trial
%   lckDI   = Discrimination Index ((in - out) / (in + out)) per trial
%
% Created 7/7/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    rzPos           % in cm, [1x2] of rZ positions
    rzLen = 30      % in cm
end

rzPos = rzPos ./ 100;
rzLen = rzLen ./ 100;

rz1 = [rzPos(1) - rzLen, rzPos(1)];
if rz1(1) < 0   % if RZ length puts the start of RZ in previous lap
    rz1(1,:) = [0, rzPos(1)];
    rz1(2,:) = [max(sess.pos(sess.lapInclude)) + (rzPos(1) - rzLen), max(sess.pos(sess.lapInclude))];
end
rz2 = [rzPos(2) - rzLen, rzPos(2)];

for i = 1:sess.nlaps
    tmpLckPos = sess.pos(sess.lckind(sess.lckind > sess.lapstt(i) & sess.lckind < sess.lapend(i)));
    
    for j = 1:size(rz1,1)
        tmpRZ1lcks(j,:) = histcounts(tmpLckPos,rz1(j,:)); % Account for anticipatory licks on previous lap
    end
    lckmap1(i,:) = tmpRZ1lcks;
    lckmap2(i,1) = histcounts(tmpLckPos,rz2);
end

if size(rz1,1) == 2
    lckmapCombine = lckmap1(:,1) + [0; lckmap1(1:end-1,2)];    % Offset by 1 lap to account for licks on prior lap, lose any licks prior to start of valTrial 1
else
    lckmapCombine = lckmap1;    % If both RZ within same lap
end

lckmap = [lckmapCombine, lckmap2];
lickDI = (lckmapCombine - lckmap2) ./ (lckmapCombine + lckmap2);

end