function [fhandle] = plot_datXdepth(root, datStruc, goodflag, muaflag, noiseflag, typeFlag)
%% Plots a metric, like Firing Rate, by depth of each unit
% Separates colormaps by unit classification (noise/mua/good)
% dot size indicates firing rate (normalized)
% dot color differentiates individual units within colormap/sorting group
%
% Inputs:
%   root        = root object. Must have root.tssync and root.tsb fields
%   datStruc    = 
%   goodflag    = binary, whether to plot good units
%   muaflag     = binary, whether to plot mua units
%   noiseflag   = binary, whether to plot noise units
%
% Outputs:
%   fhandle = handle to figure
%
% Created 8/21/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    datStruc  = 0   %struct containing metrics about units
    goodflag  = 1   % Plot good/mua/noise units
    muaflag   = 1
    noiseflag = 1
    typeFlag  = 1   % 1 = FR; 2 = Theta phase;
end

ngood   = length(root.good);
nnoise  = length(root.noise);
nmua    = length(root.mua);
spBins  = 0:30:max(root.info.depth);

try
    nshanks = numel(unique(root.info.shankID));
catch
    warning('shankID not specified, please add shankID to root.info')
    return
end

if typeFlag == 1
    % normfr      = root.info.fr / max(root.info.fr);
    maxDat      = max(root.info.fr);
    datGood     = root.info.fr(root.goodind);
    datMUA      = root.info.fr(root.muaind);
    datNoise    = root.info.fr(root.noiseind);
    xlabstr     = 'Avg Firing Rate (Hz)';
elseif typeFlag == 2
    maxDat      = 360;
    for i = 1:length(root.good)
        datGood(i) = rad2deg(datStruc.thetastats(i).ang);
        datSig(i) = datStruc.thetastats(i).p < 0.05;
    end
    datGood = datGood' + 180;
    xlabstr     = 'Theta Angle (180 = trough)';
elseif typeFlag == 3
    maxDat      = max(datStruc.trueSI);
    datGood     = datStruc.trueSI;
    datSig      = datStruc.sig <= 0.05;
    % binSig      = histcounts(datGood(datSig),spBins); %Do accumulator instead
    xlabstr     = 'Spatial Information (bits/spike)';
elseif typeFlag == 4
    maxDat      = max(datStruc.optoPkT).*1000;
    datGood     = datStruc.optoPkT'.*1000;
    xlabstr     = 'First post-opto peak (ms)';    
end

cmapgood    = sky(ngood);
cmapmua     = autumn(nmua+30);
cmapnoise   = gray(nnoise+30);

legCell = {};
legCt   = 1;

fhandle = figure; thandle = tiledlayout(1,1); 
ax1 = axes(thandle); hold on;
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])
if noiseflag
    sn = scatter(ax1, datNoise + root.info.shankID(root.noiseind).*maxDat, root.info.depth(root.noiseind), 'filled');
    sn.CData = cmapnoise(1:nnoise,:);
    sn.AlphaData = ones(length(root.noise),1)*0.5; sn.MarkerFaceAlpha = 'flat';
    legCell(legCt) = {'Noise'}; legCt = legCt +1;
end
if muaflag
    sm = scatter(ax1, datMUA + root.info.shankID(root.muaind).*maxDat, root.info.depth(root.muaind), 'filled');
    sm.CData =  cmapmua(1:nmua,:);
    sm.AlphaData = ones(length(root.mua),1)*0.5; sm.MarkerFaceAlpha = 'flat';
    legCell(legCt) = {'MUA'}; legCt = legCt +1;
end
if goodflag
    if typeFlag == 2 | typeFlag == 3
        tmpID = root.info.shankID(root.goodind);
        tmpD = root.info.depth(root.goodind);
        sgnsig = scatter(ax1, datGood(~datSig) + tmpID(~datSig).*maxDat, tmpD(~datSig), 'filled');
        sgnsig.CData = autumn(sum(~datSig));
        sgnsig.AlphaData = ones(length(root.good(~datSig)),1)*0.5; sgnsig.MarkerFaceAlpha = 'flat';
        sg = scatter(ax1, datGood(datSig) + tmpID(datSig).*maxDat, tmpD(datSig), 'filled');
        sg.CData = cmapgood(datSig',:);
        % sg.AlphaData = ones(length(root.good(datSig)),1)*0.5; sg.MarkerFaceAlpha = 'flat';
        legCell(legCt) = {'Non-Sig.'}; legCell(legCt+1) = {'Sig.'}; legCt = legCt +2;
    else
        sg = scatter(ax1, datGood + root.info.shankID(root.goodind).*maxDat, root.info.depth(root.goodind), 'filled');
        sg.CData = cmapgood;
        if legCt == 1
            sg.AlphaData = ones(length(root.good),1)*0.5; sg.MarkerFaceAlpha = 'flat';
        end
        legCell(legCt) = {'Good'}; legCt = legCt +1;
    end
end

xlim([0 maxDat*nshanks]);
if min(root.info.depth) < 500
    ylim([0 max(root.info.depth)]);
else
    ylim([min(root.info.depth) max(root.info.depth)+100]);
end
xticks(0:maxDat/3:maxDat*nshanks)
xticklabels([round(repmat([0 maxDat/3 2*maxDat/3],1,nshanks)), round(maxDat)])
ylabel('Distance from tip (um)'); 
xlabel(xlabstr)
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

% ax2 = axes(thandle); hold on;
% ax2.XAxisLocation = 'top';
% ax2.YTick = [];
% for i = 1:nshanks-1
%     plot(ax2,[i*maxDat i*maxDat], [min(root.info.depth) max(root.info.depth)], 'k--')
% end
% 
% xlim([0 maxDat*nshanks]);
% 
% ax2.Color = 'none';
% ax1.Box = 'off';
% ax2.Box = 'off';
% 
% xticks(maxDat/2:maxDat:nshanks*maxDat)
% for i = 1:nshanks
%     xtickcell{i} = sprintf('Shank %d', i-1);
% end
% xticklabels(xtickcell)

end