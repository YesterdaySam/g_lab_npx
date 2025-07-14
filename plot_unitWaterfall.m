function [fhandle,frMapSort,sortInd] = plot_unitWaterfall(frMapRaw,binedges,useSort,plotflag)
%% Plots the avg binned firing rate of all units in frMapRaw
%
% Inputs:
% frMapRaw = MxN matrix of M unit avg. firing rates in N bins
% binedges = 1xN+1 array of binedges used to creat frMapRaw
% useSort = Nx1 array of indices from previous sort to use instead of
%   default sorting based on max bin position
% plotflag = binary of whether to plot the output
%
% Outputs:
% fhandle = handle to figure
% frMapSort = sorted, normalized heatmap, MxN
% sortInd = index of sorted positions
%
% Created 7/14/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    frMapRaw        %Matrix of binned FR (Hz)
    binedges
    useSort = 0
    plotflag = 1    %binary
end

nBins = size(frMapRaw,2);
nUnits = size(frMapRaw,1);

unitMax = max(frMapRaw,[],2);
frMapNorm = frMapRaw ./ repmat(unitMax,[1, nBins]);
% lininds = find(frMapNorm == 1); %Returns linear indexing of MxN matrix
% [is,js] = ind2sub(size(frMapNorm), lininds);

if length(useSort) ~= 1
    sortInd = useSort;
else
    for i = 1:nUnits
        maxBin(i) = find(frMapNorm(i,:) == 1);
    end
    [~,sortInd] = sort(maxBin);
end

frMapSort = frMapNorm(sortInd,:);

if plotflag
    fhandle = figure; hold on; axis square
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    imagesc(frMapSort,[prctile(frMapSort,1,'all'), prctile(frMapSort,98,'all')]);
    % plot([2 2],[0 nUnits+1],'r--','LineWidth',2)
    plot([0 nBins+1],[0 nUnits+1],'k--','LineWidth',2)
    colormap("parula")
    cbar = colorbar; clim([0 0.98]);
    xticks(1:10:nBins+1)
    xticklabels(binedges(1:10:end)*100)
    xlim([0.15 nBins+0.85]); ylim([0.15 nUnits+0.85])
    xlabel('Position (cm)');
    ylabel('Unit #'); ylabel(cbar,'FR (Normalized)','FontSize',12,'Rotation',90)
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end
