function [fhandle] = plot_lckheatmap(sess,dbnsz)

arguments 
    sess
    dbnsz = 0.03
end

[binedges,~,lckrate] = plot_lickpos(sess,dbnsz,0);

fhandle = figure;
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.45])
fixRatio(fhandle);
imagesc(binedges*100,sess.valTrials,lckrate,[prctile(lckrate,1,'all'), prctile(lckrate,99,'all')]);
colormap("turbo")
cbar = colorbar; clim([0 prctile(lckrate,99,'all')]);
xlabel('Position'); % xlim([0 200])
% xticks(1:30:length(binedges)); xticklabels(binedges(1:30:length(binedges))*100);
% yticks(30:30:size(lckrate,1));
ylabel('Trial #'); ylabel(cbar,'Licks/s','FontSize',12,'Rotation',90)
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

end