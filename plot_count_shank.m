function [fhandle1, fhandle2] = plot_count_shank(root)

%% Summarize counts by shank and depth
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

end