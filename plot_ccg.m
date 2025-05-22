function [fhandle] = plot_ccg(root,unit1,unit2,bnsz,winlen,plotflag)
%% Calculates the Cross Correlogram between two units
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit1 = cluster_id
% unit2 = cluster_id
% bnsz = milliseconds bin size
% winlen = milliseconds - half length of CCG (1 winlen on either side of 0)
% plotflag = binary of whether to plot the output
%
% Outputs:
% fhandle = handle to figure
%
% Created 3/24/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    unit1
    unit2
    bnsz     = 1      % msec
    winlen   = 250    % msec
    plotflag = 1
end

spkts1 = root.ts(root.cl == unit1);
spkts2 = root.ts(root.cl == unit2);
binedges = (-winlen:bnsz:winlen)/1000;
halfInd = round(winlen/bnsz)+1;

countMat = zeros(length(spkts1),length(binedges)-1);

for i = 1:length(spkts1)
    countMat(i,:) = histcounts(spkts2,binedges+spkts1(i));
    % countMat(i,halfInd) = countMat(i,halfInd) - 1;
end

if plotflag
    fhandle = figure; hold on
    bar(1000*(binedges(1:end-1)+bnsz/1000/2),sum(countMat),'FaceColor','b');
    tmplims = ylim;
    plot([0 0],[0 tmplims(2)],'k--')
    xlabel(['Lag (ms) to unit ' num2str(unit1)])
    ylabel(['Unit ' num2str(unit2) ' Spike Count'])
    set(gca,'FontSize',12,'FontName','Arial')
end

end