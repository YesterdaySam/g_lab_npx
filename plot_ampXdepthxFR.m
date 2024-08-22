function [fhandle] = plot_ampXdepthxFR(root, goodflag, muaflag, noiseflag)
%% Plots the amplitude by depth of each unit
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

normfr      = root.info.fr / max(root.info.fr) * 500;
cmapgood    = sky(ngood);
cmapmua     = autumn(nmua+30);
cmapnoise   = gray(nnoise+30);

legCell = {};
legCt   = 1;

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])
if noiseflag
    sn = scatter(root.info.amp(root.noiseind), root.info.depth(root.noiseind), normfr(root.noiseind), cmapnoise(1:nnoise,:), 'filled');
    sn.AlphaData = ones(length(root.noise),1)*0.5; sn.MarkerFaceAlpha = 'flat';
    legCell(legCt) = {'Noise'}; legCt = legCt +1;
end
if muaflag
    sm = scatter(root.info.amp(root.muaind), root.info.depth(root.muaind), normfr(root.muaind), cmapmua(1:nmua,:), 'filled');
    sm.AlphaData = ones(length(root.mua),1)*0.5; sm.MarkerFaceAlpha = 'flat';
    legCell(legCt) = {'MUA'}; legCt = legCt +1;
end
if goodflag
    sg = scatter(root.info.amp(root.goodind), root.info.depth(root.goodind), normfr(root.goodind), cmapgood, 'filled');
    if legCt == 1
        sg.AlphaData = ones(length(root.good),1)*0.5; sg.MarkerFaceAlpha = 'flat';
    end
    legCell(legCt) = {'Good'}; legCt = legCt +1;
end

xlim([0 max(root.info.amp)]);
ylim([0 max(root.info.depth)]);
ylabel('Distance from tip (um)'); 
xlabel('Amplitude (uV)')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

end