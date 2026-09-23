function [fhandle] = plot_postProbXTime(decodeI,twin)

indStt = find(decodeI.newT > twin(1),1);
indEnd = find(decodeI.newT > twin(2),1);

fhandle = figure; hold on;
set(gcf,'Units','normalized','Position',[0 0.65 1 0.15])
imagesc(decodeI.dMat(:,indStt:indEnd))
clim([0,1])
axis xy
colormap(flipud(gray))
plot(decodeI.rPos(indStt:indEnd) / max(decodeI.rPos) * size(decodeI.dMat,1),'r')
xticks(linspace(0,indEnd-indStt,4));
xticklabels(round(linspace(decodeI.newT(indStt),decodeI.newT(indEnd),4),2));
yticks(linspace(1,size(decodeI.dMat,1),4));
yticklabels(round(linspace(0,max(decodeI.rPos),4),2));
xlabel('Time (sec)')
ylabel('Position (cm)')
ylim([0 size(decodeI.dMat,1)])
xlim([0 indEnd - indStt])
set(gca,'FontSize',16,'FontName','Arial')
end