function [] = plotPhysCompact(root,sess,sdir,overwrite,dType,rastF,velF,avgposF,trialheatmapF,templateF,perioptoF)
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
    rastF               = 1 % Include raster across trials
    velF                = 1 % Include velocity binned firing rate
    avgposF             = 1 % Include averaged spatial firing rate
    trialheatmapF       = 1 % Include heatmap across trials
    templateF           = 1 % Include template waveform 
    perioptoF           = 0 % Include peri-opto pulse spikes plot
end

cd(sdir)

if contains(dType,'good')
    if isempty(dir('ephysPlots_good')) || overwrite
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

nPlots = sum([rastF,velF,avgposF,trialheatmapF,templateF,perioptoF]);

for i = 1:nUnits

    if contains(dType,'good')
        cc = root.good(i);
    elseif contains(dType,'mua')
        cc = root.mua(i);
    end

    if sum(root.cl == cc) < 100
        disp(['Unit ' num2str(cc) ' contains ' num2str(sum(root.cl == cc)) ' spikes, skipping.'])
        continue
    end

    if isempty(dir(['unit' num2str(cc) '_summary.png'])) | overwrite == 1

        if rastF
            tmpraster = plot_trialraster(root,cc,sess);
        end

        if velF
            [~,~,~,tmpfrvel] = plot_frXvel(root,cc,sess);
            title('');
        end

        if avgposF
            [~,~,tmpfrpos] = plot_frXpos(root,cc,sess);
            title('');
        end

        if trialheatmapF
            try
                tmpheatmap = plot_trialHeatmap(root,cc,sess);
            catch
                disp(['heatmap failed for unit ' num2str(cc)])
            end
        end

        if templateF
            tmpWF = plot_templateWF(root,cc);
        end

        if perioptoF
            [~,~,tmpfropto] = plot_frXopto(root,cc,sess,0.01,0.1);
            title('');
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
        end
        tblind = find(root.info.cluster_id == cc);
        title(tcl,['unit ' num2str(cc) ' Shank ' num2str(root.info.shankID(tblind)), ' Depth ', num2str(root.info.depth(tblind)) 'um'])
        saveas(newfig, ['unit' num2str(cc) '_shank' num2str(root.info.shankID(tblind)) '_depth' num2str(root.info.depth(tblind)) '_summary'], 'png')
        close all
    end
end