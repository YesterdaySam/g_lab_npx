function [fhandle, binedges, lckmap, lckDI] = plot_lickDiscrim(sess,rzPos,rzLen,nBlock)
%% Plot Discrimiation Index of licks over trial blocks
% Inputs
%   sess    = struct from importBhvr.m
%   rzPos   = [1x2] in cm of [RZ position, comparison position]
%   rzLen   = double of reward zone length in cm
%
% Outputs
%   lckmap  = lick counts by trial
%   lckDI   = Discrimination Index ((in - out) / (in + out)) per trial
%
% Created 7/7/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    rzPos           % in cm, [1x2] of rZ positions
    rzLen  = 30     % in cm
    nBlock = 10     % nTrials per block
end

[lckmap,lckDI] = get_lickDiscrim(sess,rzPos,rzLen);

binedges = 1:nBlock:length(lckDI);

for i = 1:length(binedges)-1
    % tmp = lckDI > binedges(i) & lckDI < binedges(i+1);
    bnDI(i) = mean(lckDI(binedges(i):binedges(i+1)),'omitnan');
end

% [~,~,loc]=histcounts(1:length(lckDI),binedges);
% bnDI = accumarray(loc(:),lckDI) ./ accumarray(loc(:),1);

fhandle = figure; hold on
plot(lckDI,'k-*','LineWidth',1)
plot(binedges(1:end-1)+0.5*nBlock,bnDI,'r-o','LineWidth',2)
legend({'Trial',[num2str(nBlock) ' trial avg.']},'Location','best')
xlabel('Trial #'); ylabel('Lick DI ((RZ - AZ) / (RZ + AZ)')
set(gca,'FontSize',12,'FontName','Arial')


end