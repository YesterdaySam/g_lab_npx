function [fhandle1, fhandle2] = plot_count_shank(root, region)
%% Summarize Good unit counts by shank and depth
% Inputs
%   root = root struct with root.info.lyrID
%   region = 1 (CA1/subiculum); 2 (EC)
% Outputs
%   fhandle1 = horizontal bar for count by depth
%   fhandle2 = vertical bar for count by shank
%
% Modified 5/27/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

if region == 1
    fhandle1 = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.2 0.1 0.6])
    tmpedges = min(root.info.depth):20:max(root.info.depth);
    binCounts = histcounts(root.info.depth(root.goodind),tmpedges);
    barh(tmpedges(1:end-1)+10,binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410])
    ylim([0 max(root.info.depth)])
    xlabel('Good Counts')
    ylabel('Distance from tip (um)');
    set(gca,'FontSize',12,'FontName','Arial')

    fhandle2 = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.15])
    tmpedges = 0:1:length(unique(root.info.shankID));
    binCounts = histcounts(root.info.shankID(root.goodind),tmpedges);
    b = bar(tmpedges(1:end-1),binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410]);
    set(get(b,'Parent'),'ydir','reverse')
    ylim([0 max(binCounts)])
    xlim([-0.5 max(root.info.shankID)+0.5])
    xlabel('Shank ID')
    ylabel('Good Counts');
    set(gca,'FontSize',12,'FontName','Arial')
    text(tmpedges(1:end-1),binCounts,num2str(binCounts'),'vert','bottom','horiz','center','FontSize',12);
    % Add counts in layer by shank
    if any(ismember(root.info.Properties.VariableNames,'lyrID'))
        lyrCounts = histcounts(root.info.shankID(root.goodind & root.info.lyrID == 1),tmpedges);
        b2 = bar(tmpedges(1:end-1),lyrCounts,'FaceColor',[0.45 0.45 0.75],'EdgeColor',[0.25 0.25 0.5]);
        set(get(b2,'Parent'),'ydir','reverse')
        text(tmpedges(1:end-1),lyrCounts,num2str(lyrCounts'),'vert','bottom','horiz','center','FontSize',12);
    end

elseif region == 2
    fhandle1 = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.2 0.1 0.6])
    tmpedges = min(root.info.depth):20:max(root.info.depth);
    goodDepths = root.info.depth(root.goodind);
    binCounts = histcounts(goodDepths,tmpedges);
    lyrs = unique(root.info.lyrID(root.goodind));

    for i = 1:length(lyrs)
        lyrCounts(i,:) = histcounts(goodDepths(root.info.lyrID(root.goodind) == lyrs(i)),tmpedges);
    end
    barh(tmpedges(1:end-1)+10,binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410])
    bh1 = barh(tmpedges(1:end-1)+10,lyrCounts,'stacked');
    bh1(1).FaceColor = [1 0 1];
    bh1(2).FaceColor = [0 1 1];
    bh1(3).FaceColor = [.5 .5 .5];
    ylim([0 max(root.info.depth)])
    xlabel('Good Counts')
    ylabel('Distance from tip (um)');
    set(gca,'FontSize',12,'FontName','Arial')

    fhandle2 = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.2 0.08 0.4])
    binCounts = histcounts(root.info.lyrID(root.goodind),[2,3,5,6]);
    b = bar(0,binCounts,'stacked');
    b(1).FaceColor = [1 0 1];
    b(2).FaceColor = [0 1 1];
    b(3).FaceColor = [.5 .5 .5];
    text([0, 0, 0],cumsum(binCounts),num2str(binCounts'),'vert','top','horiz','center','FontSize',12);
    ylim([0 sum(binCounts)])
    yticks(sum(binCounts))
    xticks(unique(root.info.shankID))
    xlim([-0.5 max(root.info.shankID)+0.5])
    xlabel('Shank ID')
    ylabel('Good Counts');
    set(gca,'FontSize',12,'FontName','Arial')

end

end