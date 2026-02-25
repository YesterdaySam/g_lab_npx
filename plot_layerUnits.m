function [fhandle] = plot_layerUnits(root, goodflag, muaflag, noiseflag)
%% Plots the firing rate by depth of each unit
% Separates colormaps by unit classification (noise/mua/good)
% dot size indicates firing rate (normalized)
% dot color differentiates individual units within colormap/sorting group
%
% Inputs:
%   root        = root object. Must have root.tssync and root.tsb fields
%   goodflag    = binary, whether to plot good units
%   muaflag     = binary, whether to plot mua units
%   noiseflag   = binary, whether to plot noise units
%
% Outputs:
%   fhandle = handle to figure
%
% Created 4/3/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    goodflag  = 1   % Plot good/mua/noise units
    muaflag   = 0
    noiseflag = 0
end

ngood   = length(root.good);
nnoise  = length(root.noise);
nmua    = length(root.mua);

try
    nshanks = numel(unique(root.info.shankID));
catch
    warning('shankID not specified, please add shankID to root.info')
    return
end

normfr      = root.info.fr / max(root.info.fr);
maxfr       = max(root.info.fr);
cmapgood    = sky(ngood);
cmapmua     = autumn(nmua+30);
cmapnoise   = gray(nnoise+30);

legCell = {};
legCt   = 1;

fhandle = figure; thandle = tiledlayout(1,1); 
ax1 = axes(thandle); hold on;
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])

% if noiseflag
%     % sn = scatter(root.info.fr(root.noiseind), root.info.depth(root.noiseind), 'filled');
%     sn = scatter(ax1, root.info.fr(root.noiseind) + root.info.shankID(root.noiseind).*maxfr, root.info.depth(root.noiseind), 'filled');
%     sn.CData = cmapnoise(1:nnoise,:);
%     sn.AlphaData = ones(length(root.noise),1)*0.5; sn.MarkerFaceAlpha = 'flat';
%     legCell(legCt) = {'Noise'}; legCt = legCt +1;
% end
% if muaflag
%     % sm = scatter(root.info.fr(root.muaind), root.info.depth(root.muaind), 'filled');
%     sm = scatter(ax1, root.info.fr(root.muaind) + root.info.shankID(root.muaind).*maxfr, root.info.depth(root.muaind), 'filled');
%     sm.CData =  cmapmua(1:nmua,:);
%     sm.AlphaData = ones(length(root.mua),1)*0.5; sm.MarkerFaceAlpha = 'flat';
%     legCell(legCt) = {'MUA'}; legCt = legCt +1;
% end

% Set up depth patches
for i = 1:nshanks
    patch(ax1,[0 maxfr maxfr 0] + (i-1)*maxfr, [root.lyrbounds(1,i) root.lyrbounds(1,i) root.lyrbounds(2,i) root.lyrbounds(2,i)], [0 0 0], 'FaceAlpha', 0.05, 'EdgeColor', [1 1 1], 'HandleVisibility', 'off');
end

if goodflag
    % sg = scatter(root.info.fr(root.goodind), root.info.depth(root.goodind), 'filled');
    tmpPYLyr = root.goodind & root.info.lyrID == 1 & root.info.uType == 1;
    tmpINLyr = root.goodind & root.info.lyrID == 1 & root.info.uType == 0;
    tmpPYOut = root.goodind & root.info.lyrID ~= 1 & root.info.uType == 1;
    tmpINOut = root.goodind & root.info.lyrID ~= 1 & root.info.uType == 0;

    sgPL = scatter(ax1, root.info.fr(tmpPYLyr) + root.info.shankID(tmpPYLyr).*maxfr, root.info.depth(tmpPYLyr), 'b', 'filled');
    sgIL = scatter(ax1, root.info.fr(tmpINLyr) + root.info.shankID(tmpINLyr).*maxfr, root.info.depth(tmpINLyr), 'r', 'filled');
    sgPO = scatter(ax1, root.info.fr(tmpPYOut) + root.info.shankID(tmpPYOut).*maxfr, root.info.depth(tmpPYOut), 'filled', 'MarkerFaceColor', [.5 .5 1], 'MarkerEdgeColor', [.5 .5 1]);
    sgIO = scatter(ax1, root.info.fr(tmpINOut) + root.info.shankID(tmpINOut).*maxfr, root.info.depth(tmpINOut), 'filled', 'MarkerFaceColor', [1 .5 .5], 'MarkerEdgeColor', [1 .5 .5]);
    % sg.CData = cmapgood;
    sgPO.AlphaData = ones(sum(tmpPYOut),1)*0.5; sgPO.MarkerFaceAlpha = 'flat';
    sgIO.AlphaData = ones(sum(tmpINOut),1)*0.5; sgIO.MarkerFaceAlpha = 'flat';
    % legCell(legCt) = {'Good'}; legCt = legCt +1;
    legCell = {'Good Pyr','Good IN'};
end

xlim([0 maxfr*nshanks]);
if min(root.info.depth) < 500
    ylim([0 max(root.info.depth)]);
else
    ylim([min(root.info.depth) max(root.info.depth)+100]);
end
xticks(0:maxfr/3:maxfr*nshanks)
xticklabels([round(repmat([0 maxfr/3 2*maxfr/3],1,nshanks)), round(maxfr)])
ylabel('Distance from tip (um)'); 
xlabel('Avg Firing Rate (Hz)')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

ax2 = axes(thandle); hold on;
ax2.XAxisLocation = 'top';
ax2.YTick = [];
for i = 1:nshanks-1
    plot(ax2,[i*maxfr i*maxfr], [min(root.info.depth) max(root.info.depth)], 'k--')
end

xlim([0 maxfr*nshanks]);

ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

xticks(maxfr/2:maxfr:nshanks*maxfr)
for i = 1:nshanks
    xtickcell{i} = sprintf('Shank %d', i-1);
end
xticklabels(xtickcell)

end