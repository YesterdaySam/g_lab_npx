function [fhandle,frMapSort,sortInd] = plot_unitWaterfall(frMapRaw,binedges,useSort,plotflag,cbarF)
%% Plots the avg binned firing rate of all units in frMapRaw
%
% Inputs:
%   frMapRaw = MxN matrix of M unit avg. firing rates in N bins
%   binedges = 1xN+1 array of binedges used to creat frMapRaw
%   useSort = Nx1 array of indices from previous sort to use instead of
%       default sorting based on max bin position
%   plotflag = binary of whether to plot the output
%
% Outputs:
%   fhandle = handle to figure
%   frMapSort = sorted, normalized heatmap, MxN
%   sortInd = index of sorted positions
%
% Created 7/14/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    frMapRaw        % Matrix of binned FR (Hz)
    binedges
    useSort = 0
    plotflag = 1    % binary
    cbarF = 1
end

nBins = size(frMapRaw,2);
nUnits = size(frMapRaw,1);

unitMax = max(frMapRaw,[],2);
frMapNorm = frMapRaw ./ repmat(unitMax,[1, nBins]);

if length(useSort) ~= 1
    sortInd = useSort;
else
    for i = 1:nUnits
        tmpbns = find(frMapNorm(i,:) == 1); % In case of multiple peak normalized bins
        maxBin(i) = tmpbns(1);
    end
    [~,sortInd] = sort(maxBin);
end

frMapSort = frMapNorm(sortInd,:);

if plotflag
    fhandle = figure; hold on;
    set(gcf,'units','normalized','position',[0.4 0.35 0.25 0.5])
    set(gca,'Position',[0.13 0.11 0.75 0.815])
    imagesc(frMapSort,[prctile(frMapSort,1,'all'), prctile(frMapSort,98,'all')]);
    % plot([0 nBins+1],[0 nUnits+1],'k--','LineWidth',2)
    colormap("parula")
    if cbarF
        cbar = colorbar('Ticks',[0 1]); clim([0 1]);
        ylabel(cbar,'FR (Normalized)','FontSize',12,'Rotation',-90)
        set(cbar,'Position',[cbar.Position(1)+0.1, cbar.Position(2), 0.05, 0.25])
        ylabel('Unit');
    else
        set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.5])
        set(gca,'Position',[0.11 0.11 0.8 0.815])
    end
    xticks([1, round(nBins/2), nBins])
    xticklabels([binedges(1)*100, binedges(round(nBins/2))*100, binedges(nBins)*100])
    xlim([0.15 nBins+0.85]); ylim([0.15 nUnits+0.85])
    yticks([1,nUnits]);
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end
