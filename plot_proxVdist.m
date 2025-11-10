function [fhandle] = plot_proxVdist(dat,prox,dist)

vColors2 = [0.063 0.322 0.255; 0.353 0.612 0.223];  % Venusaur greens

nUnits = [size(dat(prox),1) size(dat(dist),1)];
bardat = [mean(dat(prox)); mean(dat(dist))];
semdat = [std(dat(prox))/sqrt(nUnits(1)); std(dat(dist))/sqrt(nUnits(2))];

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
b = bar([1.15 2.15],bardat,0.3,'FaceColor','flat','BarWidth',0.5);
b.CData = vColors2;
errorbar([1.15 2.15],bardat,semdat,'k.')
violinplot([dat(prox); dat(dist)],[1*ones(nUnits(1),1); 2*ones(nUnits(2),1)], 'ViolinColor',vColors2,'HalfViolin','left');
xlim([0.5 2.5]); ylim([0 max([dat(prox); dat(dist)],[],'all')])
xticks(1:2); xticklabels({'Proximal', 'Distal'})
set(gca,'FontSize',12,'FontName','Arial')

end
