function [] = plotEphysBasics(root,sess,sdir,overwrite,rastFlag,velFlag,avgposFlag,trialheatmapFlag)
%% Plot and save some neural analyses
%
% Inputs:
% root      = root struct from loadKS
% sdir      = string specifying path where subfolders of plots (good and mua will be saved)
% overwrite = 0 or 1 to overwrite existing plots in the dir
%
% Created 9/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root                    % root struct
    sess                    % sess struct
    sdir                    % Directory to save subplots
    overwrite           = 0 % Overwrite old plots
    rastFlag            = 0 % Include raster across trials
    velFlag             = 0 % Include velocity binned firing rate
    avgposFlag          = 0 % Include averaged spatial firing rate
    trialheatmapFlag    = 0 % Include heatmap across trials
end

cd(sdir)

if isempty(dir('*_good')) || overwrite
    mkdir('ephysPlots_good')
end
cd('ephysPlots_good')

for i = 1:length(root.good)

    cc = root.good(i);

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

    close all
end