function [binFR,binedges,zBins,zFail,sdMuNorm] = get_presence(root,unit,sess,bnsz,zthresh,plotflag)
%% Returns
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% bnsz = size of time bins, default 60sec
% zthresh = threshold of z-score below which aberrant bins will be flagged
%
% Outputs:
% binCts = spike counts per bin
% binedges = bin edges in seconds
% zBins  = z-scored FR
% zFail  = Percentage of bins falling outside zThresh
% fhandle = handle to figure
%
% Created 3/20/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    unit
    sess
    bnsz = 1   % In seconds
    zthresh = 3 % Stdevs outside mean FR
    plotflag = 0    %Binary
end

binedges = sess.ts(root.lfp_tsb(1)):bnsz:sess.ts(root.lfp_tsb(end));    % Base bins on start/end of aligned root
spkinds = root.tsb(root.cl == unit);

binCts = histcounts(sess.ts(spkinds),binedges);
binFR = binCts ./ bnsz;
muFR = mean(binFR,'omitmissing');
sdFR = std(binFR,'omitmissing');
zBins = (binFR - muFR) ./ sdFR;
superThresh = find(abs(zBins) > zthresh);
zFail = length(superThresh)/length(binedges(2:end));
sdMuNorm = sdFR/muFR;

if plotflag
    figure; hold on
    bar(binedges(1:end-1),binFR)
    plot([binedges(1) binedges(end-1)],[muFR muFR],'k--')
    if ~isempty(superThresh)
        plot(binedges(superThresh),binFR(superThresh),'r*')
    end
    xlabel('Time')
    ylabel('FR (Hz)')
    title(['Spike Presence Unit ' num2str(unit)])
    legend(['% Aberrant Bins: ' num2str(100*zFail)])
end

end