function [sess,f2rDlt] = get_QTwin(sess, qs, beltLen)
%% Get the time indices on either side of all cues, assign to sess variable
%
% Inputs:
%   sess = session struct from importBhvr
%   qs = 0 = rand cue, 1 = fixed cue; or vector of indices for other events e.g. sess.rwdind
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
    qs = 0
    beltLen = 1.85
end

% Use either the fixed or rand cue indices, or leave as is for other event indices
qflag = qs;
if qflag == 0
    qs = sess.randCueInd;
elseif qflag == 1
    qs = sess.fixedCueInd;
end

% Set spatial window length and initialize variables
wlen = beltLen / 2; % meters
trialID = [];
f2rDlt = [];
qStt = [];
qEnd = [];
qInd = [];

% Loop through each cue presentation and get the spatial window start and stop indices
for i = 2:length(qs)-1
    tmpCue = qs(i);
    qInd = [qInd; tmpCue];
    tmplap = find(sess.lapstt < tmpCue,1,'last');   % In case of cue/lap # misalignment
    trialID = [trialID; tmplap];
    qPos = sess.pos(tmpCue);
    spWinEnd = mod(qPos + wlen, beltLen);

    % Get index of end spatial bin for this cue
    if qPos + wlen > beltLen % Find ts on next lap
        tmpPosReset = find(sess.pos(tmpCue:end) < 0.01,1,'first');
        tWinEnd = find(sess.pos(tmpCue+tmpPosReset-1:end) > spWinEnd,1,'first');
        qEnd = [qEnd; tWinEnd + tmpCue + tmpPosReset - 2];
    else  % Find ts on this lap
         qEnd = [qEnd; find(sess.pos(tmpCue:end) > spWinEnd,1,'first') + tmpCue];
    end 

    % Get index of start spatial bin for this cue
    if qPos - wlen > 0
        qStt = [qStt; find(sess.pos(1:tmpCue) < spWinEnd,1,'last')]; 
    else
        tmpPosReset = find(diff(sess.pos(1:tmpCue)) < -0.3,1,'last');
        qStt = [qStt; find(sess.pos(1:tmpPosReset) < spWinEnd,1,'last') + 1];
    end

    if isscalar(qflag)
        fqPos = sess.pos(sess.fixedCueInd(tmplap)); % Get fixed cue position
        f2rDlt = [f2rDlt; fqPos - qPos, tmplap];   % Save delta between fixed cue and other cue
    end
    if length(qStt) ~= length(qEnd)
        disp('uh oh')
    end
end

% Assign to separate variables in sess depending on flag
if qflag == 0 
    sess.rqStt = qStt;
    sess.rqEnd = qEnd;
elseif qflag == 1
    sess.fqStt = qStt;
    sess.fqEnd = qEnd;
else
    sess.qStt = qStt;
    sess.qEnd = qEnd;
end
sess.qTrialID = trialID;
sess.qInd = qInd;

% Probability map of fixed cue relative to rand cue
if isscalar(qflag)
    figure; binedges = -0.925:0.025:0.925;
    set(gcf,'units','normalized','position',[0.4 0.35 0.25 0.14])
    tmpF2R = histcounts(f2rDlt,binedges);
    tmpF2R = tmpF2R ./ sum(tmpF2R);
    plot([binedges(1:length(tmpF2R))]*100,tmpF2R, 'k')
    xlabel('Distance to Rand Cue (cm)'); ylabel('P(Fixed cue)')
    set(gca,'FontSize',12,'FontName','Arial')
end

end