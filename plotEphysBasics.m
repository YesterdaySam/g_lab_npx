function [] = plotEphysBasics(root,sess,sdir,overwrite,dType,rastFlag,velFlag,avgposFlag,trialheatmapFlag,perioptoFlag)
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
    rastFlag            = 0 % Include raster across trials
    velFlag             = 0 % Include velocity binned firing rate
    avgposFlag          = 0 % Include averaged spatial firing rate
    trialheatmapFlag    = 0 % Include heatmap across trials
    perioptoFlag        = 0 % Include peri-opto pulse spikes plot
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

    if rastFlag
        if isempty(dir(['unit' num2str(cc) '_spkraster.png'])) | overwrite == 1
            tmpraster = plot_trialraster(root,cc,sess);
            saveas(tmpraster, ['unit' num2str(cc) '_spkraster'],'png')
        end
    end

    if velFlag
        if isempty(dir(['unit' num2str(cc) '_velXFR.png'])) | overwrite == 1
            [~,~,~,tmpfrvel] = plot_frXvel(root,cc,sess);
            saveas(tmpfrvel, ['unit' num2str(cc) '_velXFR'],'png')
        end
    end

    if avgposFlag
        if isempty(dir(['unit' num2str(cc) '_posXFR_5cm.png'])) | overwrite == 1
            [~,~,tmpfrpos] = plot_frXpos(root,cc,sess);
            saveas(tmpfrpos, ['unit' num2str(cc) '_posXFR_5cm'],'png')
        end
    end

    if trialheatmapFlag
        if isempty(dir(['unit' num2str(cc) '_trialFRHeatmap_5cm.png'])) | overwrite == 1
            try
                tmphtmp = plot_trialHeatmap(root,cc,sess);
                saveas(tmphtmp, ['unit' num2str(cc) '_trialFRHeatmap_5cm'],'png')
            catch
                disp(['heatmap failed for unit ' num2str(cc)])
            end
        end
    end

    if perioptoFlag
        if isempty(dir(['unit' num2str(cc) '_optoXFR.png'])) | overwrite == 1
            [~,~,tmpfrpos] = plot_frXopto(root,cc,sess);
            saveas(tmpfrpos, ['unit' num2str(cc) '_optoXFR'],'png')
        end
    end


    close all
end