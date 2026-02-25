function [fhandle] = plotColorbar(ticks,cmap)

arguments
    ticks = [0,1]
    cmap = 'parula'
end

fhandle = figure;
set(gcf,'units','normalized','position',[0.4 0.35 0.05 0.39])
ax = axes;
cbar = colorbar(ax, 'AxisLocation','in','Position',[0.05 0.05 0.6 0.9]);
ax.Visible = 'off';
colormap(cmap)
% cbar.Limits = [ticks(1), ticks(end)];
cbar.Ticks = linspace(0,1,length(ticks));
cbar.TickLabels = ticks;
set(gca,'FontSize',16,'FontName','Arial')

end