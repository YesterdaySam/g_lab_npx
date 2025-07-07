function [fhandle] = plot_prepost(root1,sess1,root2,sess2,unit,plttype)
%% Plot the firing of a unit relative to a another variable across two epochs
%
% Inputs:
% root1     = first epoch root (split with epochStruc.m)
% sess1     = first epoch sess
% root2     = second epoch root
% sess2     = second epoch sess
% unit      = cluster ID
% plttype   = Plot type. 1 = Place; 2 = Velocity; 3 = thetaMod; 4 = Reward
%
% Created 6/19/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root1
    sess1
    root2
    sess2
    unit
    plttype     % 1 = Place; 2 = Velocity; 3 = thetaMod; 4 = Reward
end

if      plttype == 1
    [~,~,fhandlePre] = plot_frXpos(root1,unit,sess1);
    [~,~,fhandlePst] = plot_frXpos(root2,unit,sess2);
    titleStr = 'Spatial Tuning ';
elseif  plttype == 2
    [~,~,~,fhandlePre] = plot_frXvel(root1,unit,sess1);
    [~,~,~,fhandlePst] = plot_frXvel(root2,unit,sess2);
    titleStr = 'Velocity Tuning ';
elseif  plttype == 3
    [~,fhandlePre] = plot_thetaMod(root1,unit,3);
    [~,fhandlePst] = plot_thetaMod(root2,unit,3);
    titleStr = 'Theta Tuning ';
elseif  plttype == 4
    [~,~,fhandlePre] = plot_frXrwdtime(root1,unit,sess1);
    [~,~,fhandlePst] = plot_frXrwdtime(root2,unit,sess2);
    titleStr = 'Reward Tuning ';
elseif  plttype == 5
    fhandlePre = plot_trialHeatmap(root1,unit,sess1);
    fhandlePst = plot_trialHeatmap(root2,unit,sess2);
    titleStr = 'Spatial Heatmap ';
end

fhandle = figure;
set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.35])

tcl = tiledlayout(fhandle, 'flow');

figlist = [fhandlePre,fhandlePst];
titlelist = {'Familiar Goal','Novel Goal'};

if plttype ~= 3; yExtrema = [min([figlist(1).Children.YLim, figlist(2).Children.YLim]), max([figlist(1).Children.YLim, figlist(2).Children.YLim])]; end

for j = 1:numel(figlist)
    figure(figlist(j))
    title(titlelist{j},'FontWeight','normal')
    ax = gca;
    if j == 2 && plttype ~= 3; ylabel(''); end
    if j == 1 && plttype == 1 | plttype == 4; ax.YAxis(2).Label.String = ''; ax.YAxis(2).TickLabels = ''; end

    if plttype ~= 3; ax.YLim = [yExtrema(1) yExtrema(2)]; end
    ax.Parent = tcl;
    ax.Layout.Tile = j;
    close(figlist(j))
end
tblind = find(root1.info.cluster_id == unit);
title(tcl,[titleStr 'Unit ' num2str(unit) ' Shank ' num2str(root1.info.shankID(tblind)), ' Depth ', num2str(root1.info.depth(tblind)) 'um'])

end