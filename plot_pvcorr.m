function [fhandle] = plot_pvcorr(dat)
% Plot PVC heatmap with colorbar scaling [0 1] of corr values
% Inputs
%   dat = correlation data in MxM spatial bins
% Outputs:
%   fhandle = handle to figure
%
% Created 4/29/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

nBins = size(dat,1);

fhandle = figure; hold on; axis square
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
imagesc(dat,[prctile(dat,0,'all'), prctile(dat,99,'all')])
colormap("turbo")
clim([0 1]);
xlim([-0 nBins+0.5]); xticks([1 nBins]); xticklabels({'0','185'}); xlabel('Position (cm) Familiar RZ')
ylim([-0 nBins+0.5]); yticks([1 nBins]); yticklabels({'0','185'}); ylabel('Position (cm) Novel RZ')
% plot([r1posInd r1posInd],[0 nBins],'w--',[0 nBins],[r2posInd r2posInd],'w--')
plot([0 nBins],[0 nBins],'w--');
plot([0 nBins/2],[nBins/2 nBins],'w--',[nBins/2 nBins],[0 nBins/2],'w--');
% cbarLims = [0 1];
% cbar = colorbar('Ticks',cbarLims,'TickLabels',cbarLims); clim(cbarLims); % ,'position', [0.9 0.12 0.05 0.3]
% ylabel(cbar,'Pearson Correlation of PV','FontSize',12,'Rotation',90)
set(gca,'FontSize',16,'FontName','Arial','YDir','normal')

end