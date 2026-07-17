parentDir = "D:\Data\Kelton\analyses\group_analyses"; 
datT = import_xldat(parentDir,"dat_include.xlsx");
aanchalDir = 'Z:\Kelton\Analyzed Data\for_Aanchal';

overwriteF = 0;

cleanInds = datT.aanchal ~= 1;
datT(cleanInds,:) = [];  %Clean excluded sessions

region = 'sub';

saveFlag = 1;
% sbase = 'subRwdShift_';
% fname = [sbase 'data8'];

dbnsz = 0.05;
histoBnsz = 5;
binedges = 0:5:185;
binpos = 0.025:dbnsz:1.825;
wlen = 150;

clear ps stats

for i = 1:height(datT)
    try
        % === Load data ===
        cd(datT.fpath{i})

        rootfile = dir("*_root.mat");

        if datT.sess_type(i) == 1
            zsubdir = fullfile(aanchalDir,'rzfam');
        elseif datT.sess_type(i) == 2
            zsubdir = fullfile(aanchalDir,'rzshift');
        else
            zsubdir = fullfile(aanchalDir,'rzrand');
        end

        cd(zsubdir)
        if exist(rootfile.name) & overwriteF ~= 1
            continue
        end
   
        cd(datT.fpath{i})
        load(rootfile.name)
        sessfile = dir("*_session.mat");
        load(sessfile.name)

        disp(['Working on ' root.name])
        nShanks = numel(unique(root.info.shankID));

        % === Assign regions for each unit ===
        regionID = zeros(height(root.info),1);
        for j = 1:nShanks
            tmpUnits = root.info.shankID == j-1;
            if datT{i,6+j}{1} == 'ca1'
                regionID(tmpUnits) = 1;
            elseif datT{i,6+j}{1} == 'sub'
                regionID(tmpUnits) = 2;
            else
                regionID(tmpUnits) = 0;
            end
        end

        % % === Identify units to include ===
        % roiUnits = zeros(length(root.good),1);
        % for j = 1:nShanks
        %     tmpUnits = root.info.shankID == j-1;
        %     if datT{i,6+j}{1} == 'ca1'
        %         if region == 'ca1'
        %             roiUnits(tmpUnits) = 1;
        %         else
        %             roiUnits(tmpUnits) = 0;
        %         end
        %     elseif datT{i,6+j}{1} == 'sub'
        %         if region == 'sub'
        %             roiUnits(tmpUnits) = 1;
        %         else
        %             roiUnits(tmpUnits) = 0;
        %         end
        %     elseif datT{i,6+j}{1} == 'bdr'
        %         roiUnits(tmpUnits) = 0;
        %     end
        % end
        % if sum(roiUnits) == 0
        %     continue
        % end

        ca1CCs = regionID == 1 & root.info.uType == true & root.info.lyrID == 1 & root.goodind;
        subCCs = regionID == 2 & root.info.uType == true & root.info.lyrID == 1 & root.goodind;

        ca1SpkMat = binarizeUnits(root,root.info.cluster_id(ca1CCs),sess);
        subSpkMat = binarizeUnits(root,root.info.cluster_id(subCCs),sess);

        save([root.name '_spkMat_ca1_pyrs'],'ca1SpkMat')
        save([root.name '_spkMat_sub_pyrs'],'subSpkMat')

        root.regionID = regionID;
        root.ca1CCs   = ca1CCs;
        root.subCCs   = subCCs;
        root.rzloc    = [datT.rzloc1(i) datT.rzloc2(i)];

        % Save to Z drive
        cd(zsubdir)
        saveRoot(root);
        saveSess(sess);
        save([root.name '_spkMat_ca1_pyrs'],'ca1SpkMat')
        save([root.name '_spkMat_sub_pyrs'],'subSpkMat')

        % Make subdirectories to save over summary plots by region
        ca1plotsDir = fullfile(zsubdir,[root.name '_ca1_plots']);
        subplotsDir = fullfile(zsubdir,[root.name '_subiculum_plots']);
        mkdir(ca1plotsDir)
        mkdir(subplotsDir)

        % Copy existing summary plots to Z drive
        cd(fullfile(datT.fpath{i},'ephysPlots_good'))
        dirlist = dir();
        dirNames = {dirlist.name};
        units = root.info.cluster_id(subCCs);
        for j = 1:length(units)
            plotID = find(contains(dirNames,['unit' num2str(units(j)) '_']));
            copyfile(dirNames{plotID}, subplotsDir);
        end
        units = root.info.cluster_id(ca1CCs);
        for j = 1:length(units)
            plotID = find(contains(dirNames,['unit' num2str(units(j)) '_']));
            copyfile(dirNames{plotID}, ca1plotsDir);
        end

    catch
        disp(['failed on ' datT.fpath{i}])
    end
end

