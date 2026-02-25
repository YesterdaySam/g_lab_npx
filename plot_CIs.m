function [] = plot_CIs(xcoords,ciup,cidn,col)
    patch([xcoords,fliplr(xcoords)],[cidn,fliplr(ciup)],col,'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
end