function [fhandle,countMat,binedges] = plot_acg(root,unit,bnsz,winlen,plotflag)
%% Calculates waveform width and FWHM and automatically assigns
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster_id
% bnsz = milliseconds bin size
% winlen = milliseconds - half length of ACG (1 winlen on either side of 0)
% plotflag = binary of whether to plot the output
%
% Outputs:
% fhandle = handle to figure
%
% Created 3/24/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    unit
    bnsz     = 1      % msec
    winlen   = 250    % msec
    plotflag = 1
end

spkts = root.ts(root.cl == unit);
binedges = (-winlen:bnsz:winlen)/1000;
halfInd = round(winlen/bnsz)+1;

countMat = zeros(length(spkts),length(binedges)-1);

for i = 1:length(spkts)
    countMat(i,:) = histcounts(spkts,binedges+spkts(i));
    countMat(i,halfInd) = countMat(i,halfInd) - 1;
end

if plotflag
    fhandle = figure;
    % bar(1000*(binedges(1:end-1)+bnsz/1000/2),sum(countMat),'FaceColor','b');
    bar(1000*(binedges(halfInd:end-1)+bnsz/1000/2),sum(countMat(:,halfInd:end)),'FaceColor','b');
    xlabel('Lag (ms)')
    ylabel('Count')
    set(gca,'FontSize',12,'FontName','Arial')
else
    fhandle = 0;
end

end