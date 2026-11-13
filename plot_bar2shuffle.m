function [fhandle] = plot_bar2shuffle(dat1,dat2,shuf1,shuf2,datcols,shufcols)

arguments
    dat1
    dat2
    shuf1
    shuf2
    datcols = [.35 .35 .35; 1 .25 .25];
    shufcols = [0.5 0.5 0.5];
end

nUnits = [size(dat1,1) size(dat2,1)];
% bardat = [mean(dat1,'omitnan'); mean(dat2,'omitnan')];
% semdat = [std(dat1)/sqrt(nUnits(1)); std(dat2)/sqrt(nUnits(2))];

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
fhandle = fixRatio(fhandle);
% b = bar([1.15 2.15],bardat,0.3,'FaceColor','flat','BarWidth',0.5);
% b.CData = cols;
% errorbar([1.15 2.15],bardat,semdat,'k.')
violinplot([dat1; dat2],[1*ones(nUnits(1),1); 2*ones(nUnits(2),1)], 'ViolinColor',datcols,'HalfViolin','left');
violinplot([shuf1; shuf2],[1*ones(nUnits(1),1); 2*ones(nUnits(2),1)], 'ViolinColor',shufcols,'HalfViolin','right');
% alldat = vertcat(dat1,dat2,shuf1,shuf2);
% [~,maxind] = max(abs(alldat),[],'all');
% ylim(sort([0 alldat(maxind)]))
xlim([0.5 2.5]); ylim([0 1]);
xticks(1:2);
xticklabels({'Strong', 'Weak'});
set(gca,'FontSize',12,'FontName','Arial')
box off

end