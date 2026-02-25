function [fhandle] = plotDeltaHisto(dat1,dat2,binedges)

deltaDatDistro = histcounts(dat2-dat1,binedges);
bnsz = binedges(2) - binedges(1);

fhandle = figure; hold on
bar(binedges(1:end-1)+bnsz/2,deltaDatDistro/sum(deltaDatDistro),'FaceColor',[0.25 0.15 1],'HandleVisibility','off');
plot(mean(dat2-dat1), max(deltaDatDistro/sum(deltaDatDistro))+0.01,'v','Color',[0.25 0.15 1])
plot([0 0], [0 max(deltaDatDistro/sum(deltaDatDistro))+0.01],'k--')
legend({['mean = ' num2str(mean(dat2-dat1))]})
ylabel('Probability')
set(gca,'FontSize',12,'FontName','Arial')
end
