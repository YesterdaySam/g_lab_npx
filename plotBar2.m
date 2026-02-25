function [fhandle] = plotBar2(dat1,dat2,vColors)

arguments
    dat1
    dat2
    % vColors = [0.5 0.5 1; 0.75 0.75 1];
    vColors = [.35 .35 .35; 1 .25 .25];
end

nUnits = [size(dat1,1) size(dat2,1)];
bardat = [mean(dat1); mean(dat2)];
semdat = [std(dat1)/sqrt(nUnits(1)); std(dat2)/sqrt(nUnits(2))];

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
b = bar([1.15 2.15],bardat,0.3,'FaceColor','flat','BarWidth',0.5);
b.CData = vColors;
errorbar([1.15 2.15],bardat,semdat,'k.')
violinplot([dat1; dat2],[1*ones(nUnits(1),1); 2*ones(nUnits(2),1)], 'ViolinColor',vColors,'HalfViolin','left');
xlim([0.5 2.5]); ylim([0 max([dat1; dat2],[],'all')])
xticks(1:2); xticklabels({'Familiar', 'Novel'})
set(gca,'FontSize',12,'FontName','Arial')
box off
end


