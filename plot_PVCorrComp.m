function [fhandle] = plot_PVCorrComp(dat1, dat2, pval)
% Plot on-diagonal vs off diagonal average PVC
% Inputs
%   dat1 = Nx1 PVC values where N = mice
%   dat2 = Nx1 PVC values in another condition
%   pval = p-value for text on graph
% Outputs
%   fhandle = handle to figure
%

nMice = length(dat1);

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.2])
plot([dat1' dat2']','-o','Color',[.5 .5 .5])
errorbar([1 2],mean([dat1' dat2'],1,'omitnan'),std([dat1' dat2'],1,'omitnan')./sqrt(nMice),'k.','LineWidth',2,'CapSize',20)
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Diagonal','Off-Diag'})
ylim([-0.25 1]); text2bar(fhandle,'Mean of PV Corr.',pval);
set(gca,'FontSize',16,'FontName','Arial')

end

