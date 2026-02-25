function [fhandle,mdlparams] = plot_anatCorr(dat,anatDist,toggle,cmap)

mdl = fitlm(anatDist, dat);
ys = predict(mdl, anatDist);
[r,p] = corrcoef(dat,ys);
r = r(2,1);
p = p(2,1);

mdlparams.r = r;
mdlparams.p = p;
mdlparams.b = mdl.Coefficients{2,1}; 
mdlparams.ypreds = ys;
mdlparams.yint = predict(mdl,0);

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])

if toggle == 1
    scatter(anatDist,dat,[],anatDist,'filled')
    colormap(cmap);
    plot([min(anatDist) max(anatDist)],predict(mdl,[min(anatDist); max(anatDist)]),'Color','k','LineWidth',2)
    xlim([0 1]); ylim([0 max([dat; dat],[],'all')])
    xlabel('% Distance through subiculum')
else
    scatter(dat,anatDist,[],anatDist,'filled')
    colormap(cmap)
    plot(predict(mdl,[min(anatDist); max(anatDist)]),[min(anatDist) max(anatDist)],'Color','k','LineWidth',2)
    ylim([min(anatDist) max(anatDist)]); xlim([0 max([dat; dat],[],'all')])
    ylabel('Distance to layer center (um)')
end

end
