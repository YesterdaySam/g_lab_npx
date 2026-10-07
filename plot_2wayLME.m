function [fhandle] = plot_2wayLME(dat1,dat2,bvgrp,cols)
% Plots data from a 2-way linear model

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.2])
fhandle = fixRatio(fhandle);
for i = 1:2
    nMice = bvgrp(i).n;
    plot([0.85*ones(nMice,1) 1.15*ones(nMice,1)]'+(i-1), [dat1(bvgrp(i).bvInd)' dat2(bvgrp(i).bvInd)']','-o','Color',cols(i,:))
    errorbar([0.85 1.15]+(i-1),mean([dat1(bvgrp(i).bvInd)' dat2(bvgrp(i).bvInd)'],1,'omitnan'),std([dat1(bvgrp(i).bvInd)' dat2(bvgrp(i).bvInd)'],1,'omitnan')./sqrt(nMice),'k.','LineWidth',2,'CapSize',20)
end
xlim([0.5 2.5]); xticks([0.85 1.15 1.85 2.15]);
ylim([-0.57 1]);
set(gca,'FontSize',16,'FontName','Arial')
end
