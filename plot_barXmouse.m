function [fhandle] = plot_barXmouse(dat,vColors)

arguments
    dat
    % vColors = [0.5 0.5 1; 0.75 0.75 1];
    vColors = [.35 .35 .35; 1 .25 .25];
end

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.1 0.27])
b = bar(mean(dat),'FaceColor','flat','HandleVisibility','off');
b.CData = vColors;
errorbar(1:2,mean(dat),std(dat)/sqrt(size(dat,1)),'k.','HandleVisibility','off')
plot(dat','k-o')
xlim([0.5 2.5])
xticks(1:2); xticklabels({'Familiar', 'Novel'})
legend({'Mouse'},'location','southeast')
set(gca,'FontSize',12,'FontName','Arial')
end