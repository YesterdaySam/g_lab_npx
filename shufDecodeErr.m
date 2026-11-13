function [shufErr,shufLoc] = shufDecodeErr(decodeI,nshufs,maxPos,plotflag)
%% Compute the shuffle error and error location
% If rPos = 0.1 and dPos = 1.8, and tracklen = 1.85, dErr = 0.06
%
% Inputs:
%   rPos     = Nx1 vector of real positions
%   dPos     = Nx1 vector of decoded positions
%   tracklen = double, in meters e.g. 1.85
%
% Outputs:
%   dErr = circular error between rPos and dPos
%
% Created 10/22/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    decodeI
    nshufs = 1000;
    maxPos = 1.85;
    plotflag = 1;
end

binedges = 0:0.05:0.95;

nts = length(decodeI.newT);
minShift = round(1 * (1/diff(decodeI.newT(1:2))));     %Minimum duration to shift spike trains (e.g. 1000ms)

for i = 1:nshufs
    rShift = minShift + randi(nts - 2*minShift,1);    % Get # of tsb to shift all spike trains, but don't over-shift
    tmpDPos = circshift(decodeI.dPos,rShift);
    tmpDErr = get_circDecodeErr(decodeI.rPos,tmpDPos,maxPos);
    shufErr(i) = mean(abs(tmpDErr));
    shufLoc(i) = get_decodeErrLoc(tmpDErr,binedges);
end

if size(shufErr,1) ~= nshufs
    shufErr = shufErr';
    shufLoc = shufLoc';
end

% Visualize
if plotflag
    figure; hold on;
    histogram(shufErr,binedges,'Normalization','probability');
    plot(repmat(mean(abs(decodeI.dErr)),2,1),[0 1],'k--');
    xlim([0 maxPos/2])
    histogram(shufLoc,binedges,'Normalization','probability');
    plot(repmat(get_decodeErrLoc(decodeI.dErr),2,1),[0 1],'r--');
end
end