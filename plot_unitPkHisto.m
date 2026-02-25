function [fhandle,frMapPk] = plot_unitPkHisto(frMapRaw,binedges,matchWfl)
%% Finds the first peak location in each row of frMapRaw and plots
%    distribution as a normalized histogram
%
% Inputs
%   frMapRaw = MxN matrix of firing rates of M units over N bins 
%   binedges = scaling of bins
%   matchWfl = Default 0, whether or not to scale plot to match waterfall
%   plot
% Outputs
%   fhandle = handle to figure
%   frMapPk = MxN binary of 1st peak in each row
%
% Created 10/31/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    frMapRaw
    binedges
    matchWfl = 0
end

nBins = size(frMapRaw,2);
nUnits = size(frMapRaw,1);

unitMax = max(frMapRaw,[],2);
frMapPk = zeros(size(frMapRaw));

for i = 1:nUnits
    tmpbns = find(frMapRaw(i,:) == unitMax(i)); % In case of multiple peak normalized bins
    frMapPk(i,tmpbns(1)) = 1;
end

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.25 0.14])
plot(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(frMapPk)./sum(frMapPk,'all'),'Color',[0.25 0.15 1])
% bar(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(frMapPk)./sum(frMapPk,'all'),'FaceColor',[0.25 0.15 1])
xlim([0 max(binedges)])

if matchWfl
    set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.14])
    set(gca,'Position',[0.11 0.17 0.8 0.80])
    nEdges = nBins+1;
    xticks([binedges(1), binedges(round(nEdges/2)), binedges(nEdges)])
else
    ylabel('Probability')
end
set(gca,'FontSize',12,'FontName','Arial')
end