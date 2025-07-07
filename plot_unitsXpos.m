function [fhandle,frMapSort,sortInd] = plot_unitsXpos(root,sess,units,useSort,dbnsz,vthresh,plotflag)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% sess = session struct from importBhvr
% units = Nx1 array of cluster IDs
% useSort = Nx1 array of indices from previous sort to use instead of
%   default sorting based on max bin position
% dbnsz = size of spatial bins, default 0.04 = 4cm
% vbnsz = size of velocity bins, default 0.02m/s = 2cm/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% 
%
% Created 9/19/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    sess            %session struct
    units {double}  %Cluster IDs of the root object
    useSort = 0
    dbnsz = 0.05    %m
    vthresh = 0.04  %m/s; velocity threshold for spikes
    plotflag = 1    %binary
end

% Only use valid trials
sess2 = sess;
sess2.lapstt = sess.lapstt(sess.valTrials);
sess2.lapend = sess.lapend(sess.valTrials);

nUnits = length(units);
maxpos = max(sess2.pos(sess2.lapstt(1):sess2.lapend(1)));
binedges = 0:dbnsz:maxpos;    % Base max binsize on first valid trial
nBins = length(binedges)-1;
frMapRaw = zeros(nUnits, nBins);

for i = 1:nUnits
    % [~,tmpbnfr] = plot_frXpos(root,units(i),sess,dbnsz,vthresh,0);
    % [~,~,~,~,tmpbnfr] = get_PF(root,units(i),sess,dbnsz,vthresh);
    [~,~,~,~,~,~,tmpbnfr] = get_SI(root,units(i),sess,dbnsz,vthresh);
    % frMapRaw(i,:) = mean(tmpbnfr,1,'omitnan');
    frMapRaw(i,:) = tmpbnfr;
end

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

end
