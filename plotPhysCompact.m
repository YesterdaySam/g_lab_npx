function [] = plotPhysCompact(root,sess,sdir,overwrite,dType,rastFlag,velFlag,avgposFlag,trialheatmapFlag,perioptoFlag)
%% Plot and save some neural analyses
%
% Inputs:
% root      = root struct from loadKS
% sess      = session struct from importBhvr
% sdir      = string specifying path where subfolders of plots (good and mua will be saved)
% overwrite = 0 or 1 to overwrite existing plots in the dir
% dType     = string specifying 'good' or 'mua' for which units to plot
%
% Created 9/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root                    % root struct
    sess                    % sess struct
    sdir                    % Directory to save subplots
    overwrite           = 0 % Overwrite old plots
    dType               = 'good' % Type of units to work with
    rastFlag            = 1 % Include raster across trials
    velFlag             = 1 % Include velocity binned firing rate
    avgposFlag          = 1 % Include averaged spatial firing rate
    trialheatmapFlag    = 1 % Include heatmap across trials
    perioptoFlag        = 1 % Include peri-opto pulse spikes plot
end

cd(sdir)

if contains(dType,'good')
    if isempty(dir('*_good')) || overwrite
        mkdir('ephysPlots_good')
    end
    cd('ephysPlots_good')
    nUnits = length(root.good);
elseif contains(dType,'mua')
    if isempty(dir('*_mua')) || overwrite
        mkdir('ephysPlots_mua')
    end
    cd('ephysPlots_mua')
    nUnits = length(root.mua);
end

nPlots = sum([rastFlag,velFlag,avgposFlag,trialheatmapFlag,perioptoFlag]);

for i = 1:nUnits

    if contains(dType,'good')
        cc = root.good(i);
    elseif contains(dType,'mua')
        cc = root.mua(i);
    end

    if isempty(dir(['unit' num2str(cc) '_summary.png'])) | overwrite == 1

        if rastFlag
            tmpraster = plot_trialraster(root,cc,sess);
        end

        if velFlag
            [~,~,~,tmpfrvel] = plot_frXvel(root,cc,sess);
            title('');
        end

        if avgposFlag
            [~,~,tmpfrpos] = plot_frXpos(root,cc,sess);
            title('');
        end

        if trialheatmapFlag
            try
                tmpheatmap = plot_trialHeatmap(root,cc,sess);
            catch
                disp(['heatmap failed for unit ' num2str(cc)])
            end
        end

        if perioptoFlag
            [~,~,tmpfropto] = plot_frXopto(root,cc,sess);
            title('');
        end

        figlist = get(groot, 'Children');
        newfig = figure;
        if nPlots < 5
            set(gcf,'units','normalized','position',[0.2 0.1 0.4 0.7])
        else
            set(gcf,'units','normalized','position',[0.2 0.1 0.7 0.7])
        end

        tcl = tiledlayout(newfig, 'flow');

        for j = 1:numel(figlist)
            figure(figlist(j))
            ax = gca;
            ax.Parent = tcl;
            ax.Layout.Tile = j;
        end

        saveas(newfig, ['unit' num2str(cc) '_summary'], 'png')
        close all
    end
end