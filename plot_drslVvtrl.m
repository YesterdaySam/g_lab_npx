function [fhandle] = plot_drslVvtrl(dat,drsl,vtrl)

% vColors2 = [1 0.482 0.451; 0.87 0.255 0.255]; % Venusaur reds
vColors2 = [0.353 0.835 0.772; 0.063 0.482 0.416]; % Venusaur blues

nUnits = [size(dat(drsl),1) size(dat(vtrl),1)];
bardat = [mean(dat(drsl)); mean(dat(vtrl))];
semdat = [std(dat(drsl))/sqrt(nUnits(1)); std(dat(vtrl))/sqrt(nUnits(2))];

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
b = barh([1.15 2.15],bardat,0.3,'FaceColor','flat','BarWidth',0.5);
b.CData = vColors2;
errorbar(bardat,[1.15 2.15],semdat,'k.','horizontal')
violinplot([dat(drsl); dat(vtrl)],[1*ones(nUnits(1),1); 2*ones(nUnits(2),1)], 'ViolinColor',vColors2,'HalfViolin','left','Orientation','horizontal');
ylim([0.5 2.5]); xlim([0 max([dat(drsl); dat(vtrl)],[],'all')])
yticks(1:2); yticklabels({'Dorsal', 'Ventral'}); ytickangle(90)
set(gca,'FontSize',12,'FontName','Arial','YDir','reverse')

end
