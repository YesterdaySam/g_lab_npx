function [fhandle] = plotDeltaHisto2(datA1,datA2,datB1,datB2,ps,binedges,vColors2)

deltaDistA = histcounts(datA2-datA1,binedges);
deltaDistB = histcounts(datB2-datB1,binedges);
bnsz = binedges(2) - binedges(1);
bnctrs = binedges(1:end-1)+bnsz/2;

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])

bA = bar(bnctrs,deltaDistA/sum(deltaDistA),'FaceColor',vColors2(1,:),'FaceAlpha',0.5,'HandleVisibility','off');
bB = bar(bnctrs,deltaDistB/sum(deltaDistB),'FaceColor',vColors2(2,:),'FaceAlpha',0.5,'HandleVisibility','off');
plot(mean(datA2-datA1), max(deltaDistA/sum(deltaDistA))+0.01,'v','Color',vColors2(1,:))
plot(mean(datB2-datB1), max(deltaDistB/sum(deltaDistB))+0.01,'v','Color',vColors2(2,:))
ylims = ylim; xlims = xlim;
plot([0 0],ylims,'k--')
% plot([0 0], [0 max(deltaDistA/sum(deltaDistA))+0.01],'k--')

% text(xlims(2) - 0.5*diff(xlims), ylims(2)-.1*diff(ylims), ['mean = ' num2str(round(mean(datA2-datA1),1))], 'FontSize', 12,'Color',vColors2(1,:))
text(xlims(2) - 0.4*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(round(ps(1),3))], 'FontSize', 12,'Color',vColors2(1,:))
% text(xlims(2) - 0.5*diff(xlims), ylims(2)-.2*diff(ylims), ['mean = ' num2str(round(mean(datB2-datB1),1))], 'FontSize', 12,'Color',vColors2(2,:))
text(xlims(2) - 0.4*diff(xlims), ylims(2)-.25*diff(ylims), ['p = ' num2str(round(ps(2),3))], 'FontSize', 12,'Color',vColors2(2,:))

% legend({['mean = ' num2str(round(mean(datA2-datA1),1))], ['mean = ' num2str(round(mean(datB2-datB1),1))]})
ylabel('Probability')
set(gca,'FontSize',16,'FontName','Arial')
end
