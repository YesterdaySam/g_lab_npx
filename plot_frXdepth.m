function [fhandle] = plot_frXdepth(root, goodflag, muaflag, noiseflag)
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
% Created 8/21/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    goodflag  = 1   % Plot good/mua/noise units
    muaflag   = 1
    noiseflag = 1
end

ngood   = length(root.good);
nnoise  = length(root.noise);
nmua    = length(root.mua);

cmapgood    = sky(ngood);
cmapmua     = autumn(nmua+30);
cmapnoise   = gray(nnoise+30);

legCell = {};
legCt   = 1;

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])
if noiseflag
    sn = scatter(root.info.fr(root.noiseind), root.info.depth(root.noiseind), 'filled');
    sn.CData = cmapnoise(1:nnoise,:);
    sn.AlphaData = ones(length(root.noise),1)*0.5; sn.MarkerFaceAlpha = 'flat';
    legCell(legCt) = {'Noise'}; legCt = legCt +1;
end
if muaflag
    sm = scatter(root.info.fr(root.muaind), root.info.depth(root.muaind), 'filled');
    sm.CData =  cmapmua(1:nmua,:);
    sm.AlphaData = ones(length(root.mua),1)*0.5; sm.MarkerFaceAlpha = 'flat';
    legCell(legCt) = {'MUA'}; legCt = legCt +1;
end
if goodflag
    sg = scatter(root.info.fr(root.goodind), root.info.depth(root.goodind), 'filled');
    sg.CData = cmapgood;
    if legCt == 1
        sg.AlphaData = ones(length(root.good),1)*0.5; sg.MarkerFaceAlpha = 'flat';
    end
    legCell(legCt) = {'Good'}; legCt = legCt +1;
end

xlim([0 max(root.info.fr)]);
ylim([0 max(root.info.depth)]);
ylabel('Distance from tip (um)'); 
xlabel('Firing Rate (Hz)')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

end