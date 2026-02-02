function [sess,f2rDlt] = get_QTwin(sess, qFlag, beltLen)
%% Get the time indices on either side of all cues, assign to sess variable
%
% Inputs:
%   sess = session struct from importBhvr
%   qFlag = 0 for rand cue and 
%   beltLen = length of treadmill belt, default 1.85m
%
% Outputs:
%   sess = updated sess file with fields for the cue start and end
%       rqStt/fqStt = starts; depending on random or fixed cue type
%       rqEnd/fqEnd = ends; depending on random or fixed cue type
%
% Created 1/29/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    qFlag = 0
    beltLen = 1.85
end

% Use either the fixed or rand cue indices
if qFlag == 0
    qs = sess.randCueInd;
else
    qs = sess.fixedCueInd;
end

% Set spatial window length and initialize variables
wlen = beltLen / 2; % meters
qStt = [];
qEnd = [];
f2rDlt = [];

% Loop through each cue presentation and get the spatial window start and stop indices
for i = 2:length(qs)-1
    tmpCue = qs(i);
    tmplap = find(sess.lapstt < tmpCue,1,'last');   % In case of cue/lap # misalignment
    qPos = sess.pos(tmpCue);
    fqPos = sess.pos(sess.fixedCueInd(tmplap)); % Get fixed cue position
    spWinEnd = mod(qPos + wlen, beltLen);

    % Get index of end spatial bin for this cue
    if qPos + wlen > beltLen % Find ts on next lap
        tmpPosReset = find(sess.pos(tmpCue:end) < 0.01,1,'first');
        tWinEnd = find(sess.pos(tmpCue+tmpPosReset-1:end) > spWinEnd,1,'first');
        qEnd(i) = tWinEnd + tmpCue + tmpPosReset - 2;
    else  % Find ts on this lap
        qEnd(i) = find(sess.pos(tmpCue:end) > spWinEnd,1,'first') + tmpCue;
    end 

    % Get index of start spatial bin for this cue
    if qPos - wlen > 0
        qStt(i) = find(sess.pos(1:tmpCue) < spWinEnd,1,'last'); 
    else
        tmpPosReset = find(diff(sess.pos(1:tmpCue)) < -0.3,1,'last');
        qStt(i) = find(sess.pos(1:tmpPosReset) < spWinEnd,1,'last') + 1;
    end

    f2rDlt = [f2rDlt; fqPos - qPos, tmplap];   % Save delta between fixed cue and other cue
end

% Assign to separate variables in sess depending on flag
if qFlag == 0 
    sess.rqStt = qStt;
    sess.rqEnd = qEnd;
else
    sess.fqStt = qStt;
    sess.fqEnd = qEnd;
end

% Probability map of fixed cue relative to rand cue
figure; binedges = -0.925:0.025:0.925;
set(gcf,'units','normalized','position',[0.4 0.35 0.25 0.14])
tmpF2R = histcounts(f2rDlt,binedges);
tmpF2R = tmpF2R ./ sum(tmpF2R);
plot([binedges(1:length(tmpF2R))]*100,tmpF2R, 'k')
xlabel('Distance to Rand Cue (cm)'); ylabel('P(Fixed cue)')
set(gca,'FontSize',12,'FontName','Arial')

end