function [ripSpkPr, binedges, fhandle] = plot_frXripple(root,unit,wlen,bnsz,plotflag)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% wlen = size of window on either side of ripple peak, msec
% bnsz = size of binning for spike probability line, msec
% plotflag = binary of whether to plot the output
%
% Outputs:
% ripSpkPr = probability over entire +/- wlen of spikes falling in each bin
% binedges = temporal bin edges
% fhandle = handle to figure
%
% Created 1/10/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    wlen     = 125  %msec
    bnsz     = 5    %msec
    plotflag = 1    %binary
end

nRips = size(root.ripples,1);
wdw = wlen/1000*root.fs_lfp;
ripFRMap = zeros(nRips, wlen*2/bnsz);
spkinds = root.tsb(root.cl == unit);
% binedges = -wlen:bnsz:wlen;
binedges = -wdw:bnsz/1000*root.fs_lfp:wdw;
ripRast = [];

for i = 1:nRips
    % sigInds = [find(root.lfp_tsb == root.ripples(i,2) - wdw) find(root.lfp_tsb == root.ripples(i,2) + wdw)];
    tmpspks         = spkinds(spkinds > root.ripples(i,2) - wdw & spkinds < root.ripples(i,2) + wdw) - root.ripples(i,2);
    ripFRMap(i,:)   = histcounts(tmpspks, binedges);
    ripRast = [ripRast; tmpspks, i*ones(numel(tmpspks),1)];
end

ripSpkPr  = sum(ripFRMap) ./ sum(ripFRMap,'all');

if plotflag
    fhandle = figure; hold on;
    % set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    plot(ripRast(:,1)*1000/root.fs_lfp,ripRast(:,2),'k|')
    xlabel('Time to ripple peak (msec)')
    ylabel('Ripple #')
    currlims = ylim;
    ylim([currlims(1) currlims(2)*1.10])
    xlim([-wlen wlen])
    yyaxis right
    plot(-wlen:bnsz:wlen-1, ripSpkPr, 'r')
    ylabel('Spike Probability'); ylim([0 1])
    ax = gca;
    ax.YAxis(2).Color = 'r';
    % legend({'Spike Raster','Spike Probability'})
    set(gca,'FontSize',12,'FontName','Arial')
    tmpind = find(root.info.cluster_id == unit);
    title(['Unit ' num2str(unit) ', Shank ' num2str(root.info.shankID(tmpind)), ', Depth ', num2str(root.info.depth(tmpind)) 'um'])
end

end