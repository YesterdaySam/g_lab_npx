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
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
% vColors2 = [0.5 0.5 1; 0.75 0.75 1];
vColors2 = [.35 .35 .35; 1 .25 .25];

clear ps stats

for i = 1:height(datT)
    try
        % === Load data ===
        cd(datT.fpath{i})

        rootfile = dir("*_root.mat");

        cd(aanchalDir)
        if exist(rootfile.name) & overwriteF ~= 1
            continue
        end
   
        cd(datT.fpath{i})
        load(rootfile.name)
        sessfile = dir("*_session.mat");
        load(sessfile.name)

        disp(['Working on ' root.name])
        nShanks = numel(unique(root.info.shankID));

        % === Identify units to include ===
        roiUnits = zeros(length(root.good),1);
        for j = 1:nShanks
            tmpUnits = root.info.shankID == j-1;
            if datT{i,6+j}{1} == 'ca1'
                if region == 'ca1'
                    roiUnits(tmpUnits) = 1;
                else
                    roiUnits(tmpUnits) = 0;
                end
            elseif datT{i,6+j}{1} == 'sub'
                if region == 'sub'
                    roiUnits(tmpUnits) = 1;
                else
                    roiUnits(tmpUnits) = 0;
                end
            elseif datT{i,6+j}{1} == 'bdr'
                roiUnits(tmpUnits) = 0;
            end
        end
        if sum(roiUnits) == 0
            continue
        end

        useCCs = roiUnits & root.info.uType == true & root.info.lyrID == 1 & root.goodind;

        subSpkMat = binarizeUnits(root,root.info.cluster_id(useCCs),sess);

        save([root.name '_spkMat_sub_pyrs'],'subSpkMat')

        root.roiUnits = roiUnits;
        root.useCCs = useCCs;

        % Save to Z drive
        cd(aanchalDir)
        saveRoot(root);
        saveSess(sess);
        save([root.name '_spkMat_sub_pyrs'],'subSpkMat')
        subplotsDir = fullfile(aanchalDir,[root.name '_plots']);
        mkdir(subplotsDir)

        % Copy existing summary plots to Z drive
        cd(fullfile(datT.fpath{i},'ephysPlots_good'))
        dirlist = dir();
        dirNames = {dirlist.name};
        units = root.info.cluster_id(useCCs);
        for j = 1:length(units)
            plotID = find(contains(dirNames,['unit' num2str(units(j)) '_']));
            copyfile(dirNames{plotID}, subplotsDir);
        end
    catch
        disp(['failed on ' datT.fpath{i}])
    end
end

