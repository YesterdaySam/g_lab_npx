function [lckRate, lickDI] = get_lickDiscrim(sess,rzPos,rzLen)
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

ct = 1;
for i = sess.valTrials
    smoothlap = smooth(sess.pos(sess.lapstt(i):sess.lapend(i)),20);
    tmpLckPos = sess.pos(sess.lckind(sess.lckind > sess.lapstt(i) & sess.lckind < sess.lapend(i)));

    rz1 = [rzPos(1) - rzLen, rzPos(1)];
    if rz1(1) < 0   % if RZ length puts the start of RZ in previous lap
        rz1(1,:) = [0, rzPos(1)];
        rz1(2,:) = [max(smoothlap) + (rzPos(1) - rzLen), max(smoothlap)];
    end
    rz2 = [rzPos(2) - rzLen, rzPos(2)];
    if rz2(1) < 0   % if RZ length puts the start of RZ in previous lap
        rz2(1,:) = [0, rzPos(2)];
        rz2(2,:) = [max(smoothlap) + (rzPos(2) - rzLen), max(smoothlap)];
    end

    for j = 1:size(rz1,1)
        tmpRZ1lcks(j,:) = histcounts(tmpLckPos,rz1(j,:)); % Account for anticipatory licks on previous lap
        rz1Onn(j) = find(smoothlap > rz1(j,1),1);
        rz1Off(j) = find(smoothlap < rz1(j,2),1,'last');
        rz1Occ(j) = (rz1Off(j) - rz1Onn(j)) ./ sess.samprate;   % Get occupancy in RZ1 prep zone
    end

    for j = 1:size(rz2,1)
        tmpRZ2lcks(j,:) = histcounts(tmpLckPos,rz2(j,:)); % Account for anticipatory licks on previous lap
        rz2Onn(j) = find(smoothlap > rz2(j,1),1);
        rz2Off(j) = find(smoothlap < rz2(j,2),1,'last');
        rz2Occ(j) = (rz2Off(j) - rz2Onn(j)) ./ sess.samprate;   % Get occupancy in RZ1 prep zone
    end

    rwdInd = find(sess.rwdTrials == i,1);
    if ~isempty(rwdInd) % Add 1 lick for first lick in RZ 1 or 2, depending on reward site
        try
        rwdPos = sess.pos(sess.rwdind(rwdInd)); %Find where reward is delivered
        if rwdPos > rzPos(1) & rwdPos < rzPos(1) + 0.3 % Hardcoded 30cm RZ
            tmpRZ1lcks(1,1) = tmpRZ1lcks(1,1) + 1;
        elseif rwdPos > rzPos(2) & rwdPos < rzPos(2) + 0.3 % Hardcoded 30cm RZ
            tmpRZ2lcks(1,1) = tmpRZ2lcks(1,1) + 1;
        end
        catch
            disp('uh on')
        end
    end
    
    lckmap1(ct,:) = tmpRZ1lcks;
    lckmap2(ct,:) = tmpRZ2lcks;

    occmap1(ct,:) = rz1Occ;
    occmap2(ct,:) = rz2Occ;

    ct = ct+1;
end

if size(rz1,1) == 2
    lckmap1Combine = lckmap1(:,1) + [0; lckmap1(1:end-1,2)];    % Offset by 1 lap to account for licks on prior lap, lose any licks prior to start of valTrial 1
    occmap1Combine = occmap1(:,1) + [0; occmap1(1:end-1,2)];
else
    lckmap1Combine = lckmap1;    % If both RZ within same lap
    occmap1Combine = occmap1;
end
if size(rz2,1) == 2
    lckmap2Combine = lckmap2(:,1) + [0; lckmap2(1:end-1,2)];    % Offset by 1 lap to account for licks on prior lap, lose any licks prior to start of valTrial 1
    occmap2Combine = occmap2(:,1) + [0; occmap2(1:end-1,2)];
else
    lckmap2Combine = lckmap2;    % If both RZ within same lap
    occmap2Combine = occmap2;
end


lckmap = [lckmap1Combine, lckmap2Combine];
lckRate = lckmap ./ [occmap1Combine, occmap2Combine];
lickDI = (lckRate(:,1) - lckRate(:,2)) ./ (lckRate(:,1) + lckRate(:,2));
nanDI  = isnan(lickDI);
lickDI(nanDI) = 0;
end