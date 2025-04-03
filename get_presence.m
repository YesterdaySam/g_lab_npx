function [binCts,zFail,sdMuNorm] = get_presence(root,unit,sess,bnsz,zthresh,plotflag)
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
% zFail  = Percentage of bins falling outside zThresh
% fhandle = handle to figure
%
% Created 3/20/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    unit
    sess
    bnsz = 60   % In seconds
    zthresh = 2 % Stdevs above mean FR
    plotflag = 0    %Binary
end

binedges = sess.ts(root.lfp_tsb(1)):bnsz:sess.ts(root.lfp_tsb(end));    % Base bins on start/end of aligned root
spkinds = root.tsb(root.cl == unit);

binCts = histcounts(sess.ts(spkinds),binedges);
muCt = mean(binCts,'omitmissing');
sdCt = std(binCts,'omitmissing');
zBins = (binCts - muCt) ./ sdCt;
superThresh = find(abs(zBins) > zthresh);
zFail = length(superThresh)/length(binedges(2:end));
sdMuNorm = sdCt/muCt;

if plotflag
    figure; hold on
    bar(binedges(1:end-1)/bnsz,binCts)
    plot([binedges(1)/bnsz binedges(end-1)/bnsz],[muCt muCt],'k--')
    if ~isempty(superThresh)
        plot(binedges(superThresh)/bnsz,binCts(superThresh),'r*')
    end
    xlabel('Time (min)')
    ylabel('Spk Count')
    title(['Spike Presence Unit ' num2str(unit)])
    legend(['% Aberrant Bins: ' num2str(100*zFail)])
end

end