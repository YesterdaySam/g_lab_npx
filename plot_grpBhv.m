function [fhandle] = plot_grpBhv(dat,grp,vcolors)


ngrps = length(grp);
nCons = size(dat,2);

for i = 1:ngrps
    subdat(i).dat = dat(grp(i).bvInd,:);
    grpMean(i,:) = mean(subdat(i).dat);
    grpSEM(i,:) = std(subdat(i).dat)./sqrt(length(subdat(i).dat));
end

xcoords = 0.85:1:ngrps;
xshift = 0.3/(ngrps-1);

fhandle = figure; hold on
% set(gcf,'Units','normalized','Position',[0.1 0.45 0.376 0.186])
set(gcf,'Units','normalized','Position',[0.1 0.4 0.3333 0.1786])
% fhandle = fixRatio(fhandle);
for i = 1:ngrps
    errorbar(xcoords + xshift*(i-1), grpMean(i,:), grpSEM(i,:),'k.');
    plot(ones(size(subdat(i).dat)) .* (xcoords + xshift*(i-1)), subdat(i).dat,'o','Color', vcolors(i,:))
end
xticks([1 ngrps]); xlim([0.5 nCons+0.5]); xlabel('Condition')
set(gca,'FontSize',16,'FontName','Arial')

end