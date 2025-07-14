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

axisFlag = 0;

if      plttype == 1
    [~,~,fhandlePre] = plot_frXpos(root1,unit,sess1);
    [~,~,fhandlePst] = plot_frXpos(root2,unit,sess2);
    titleStr = 'Spatial Tuning ';
    axisFlag = 1;
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
    axisFlag = 1;
elseif  plttype == 5
    fhandlePre = plot_trialHeatmap(root1,unit,sess1);
    fhandlePst = plot_trialHeatmap(root2,unit,sess2);
    titleStr = 'Spatial Heatmap ';
elseif  plttype == 6
    [~,bnVelPre] = plot_trialvel(sess1, 0.01, 0);
    [~,bnVelPst] = plot_trialvel(sess2, 0.01, 0);
    fhandlePre = plot_vel_lck(sess1, bnVelPre, 0.01, 0.03);
    fhandlePst = plot_vel_lck(sess2, bnVelPst, 0.01, 0.03);
    titleStr = 'Velocity & Lick Rate ';
    axisFlag = 1;
elseif  plttype == 7
    [~,~,fhandlePre] = plot_frXripple(root1,unit,sess1, root1.ripRef,150);
    [~,~,fhandlePst] = plot_frXripple(root2,unit,sess2, root2.ripRef,150);
    titleStr = 'SPW-R Tuning ';
    axisFlag = 1;
end

fhandle = figure;
set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.35])

tcl = tiledlayout(fhandle, 'flow');

figlist = [fhandlePre,fhandlePst];
titlelist = {'Familiar Goal','Novel Goal'};

if plttype ~= 3
    for i = 1:length(figlist(1).Children.YAxis)
        yExtrema(i,:) = [min([figlist(1).Children.YAxis(i).Limits, figlist(2).Children.YAxis(i).Limits]),...
            max([figlist(1).Children.YAxis(i).Limits figlist(2).Children.YAxis(i).Limits])];
    end
end

for j = 1:numel(figlist)
    figure(figlist(j))
    title(titlelist{j},'FontWeight','normal')
    ax = gca;
    if j == 2 && plttype ~= 3; ax.YAxis(1).Label.String = ''; end
    if j == 1 && axisFlag; ax.YAxis(2).Label.String = ''; ax.YAxis(2).TickLabels = ''; end

    if plttype ~= 3; 
        for i = 1:length(ax.YAxis)
            ax.YAxis(i).Limits = yExtrema(i,:);
        end
    end
    ax.Parent = tcl;
    ax.Layout.Tile = j;
    close(figlist(j))
end
tblind = find(root1.info.cluster_id == unit);
if plttype ~= 6; title(tcl,[titleStr 'Unit ' num2str(unit) ' Shank ' num2str(root1.info.shankID(tblind)), ' Depth ', num2str(root1.info.depth(tblind)) 'um']); end

end