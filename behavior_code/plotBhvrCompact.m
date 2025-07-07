function [] = plotBhvrCompact(sess,sdir,overwrite,velheatmapFlag,lickFlag,rwdlickFlag,velXlickFlag,trialVelXcorrFlag)
%% Plot and save some neural analyses
%
% Inputs:
% sess      = session struct from importBhvr
% sdir      = string specifying path where subfolders of plots (good and mua will be saved)
% overwrite = 0 or 1 to overwrite existing plots in the dir
%
% Created 12/2/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess                    % sess struct
    sdir                    % Directory to save subplots
    overwrite           = 0 % Overwrite old plots
    velheatmapFlag      = 1 % Include heatmap of velocity by trial
    lickFlag            = 1 % Include lick raster by trial
    rwdlickFlag         = 1 % Include reward aligned licking
    velXlickFlag        = 1 % Include average velocity and lick rate plot
    trialVelXcorrFlag   = 1 % Include 2x heatmaps of spatially binned velocity by trial
end

cd(sdir)

nPlots = sum([velheatmapFlag,lickFlag,rwdlickFlag,velXlickFlag,trialVelXcorrFlag*2]);

if isempty(dir([sess.name(1:14) '_behavior_summary.png'])) | overwrite == 1

    if trialVelXcorrFlag
        [fig_XcorrLap, fig_XcorrTrial] = plot_lap_velCCorr(sess);
    end

    if rwdlickFlag
        try
            [fig_lickpos, fig_licktrialavg, tmpedges2, ~, tmpbnlck] = plot_lickpos(sess);
            close(fig_licktrialavg)
            title("Licks by Trial")
        catch
            disp('Reward-aligned lick raster failed')
        end
    end

    if lickFlag
        [fig_lickraster, fig_lickavg] = plot_lickraster_rwd(sess, 5, 0.1);
        close(fig_lickavg)
        title("Reward Aligned Licks")
    end

    if velheatmapFlag
        [~,tmpbnvel,fig_trialvel,fig_velavg] = plot_trialvel(sess);  % default 0.01m binsize
        close(fig_velavg)
        title("Velocity by Trial")
    end

    if velXlickFlag
        fig_vel_lick = plot_vel_lck(sess, tmpbnvel, 0.01, 0.03);
        title("Average Lickrate and Velocity")
    end

    figlist = get(groot, 'Children');
    newfig = figure;
    if nPlots < 5
        set(gcf,'units','normalized','position',[0.1 0.1 0.4 0.7])
    else
        set(gcf,'units','normalized','position',[0.1 0.1 0.7 0.7])
    end

    tcl = tiledlayout(newfig, 'flow');

    for j = 1:numel(figlist)
        figure(figlist(j))
        ax = gca;
        ax.Parent = tcl;
        ax.Layout.Tile = j;
        if ax.Title.String == "Velocity by Trial"
            colormap(ax,'sky')
        end
    end
    title(tcl,replace(sess.name(1:end-8),'_',' '))
    saveas(newfig, [sess.name(1:14) '_behavior_summary.png'], 'png')
    close all
end
end