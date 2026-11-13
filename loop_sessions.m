function [output] = loop_sessions(datT)

for i = 1:height(datT)

    % === Load data ===
    cd(datT.fpath{i})
    [root,sess] = get_rootSessDat(datT.fpath{i},1,1,0);
    epochfile = dir("*_dat.mat");
    load(epochfile.name)

    % === Plot Phys Compact - Good, in-layer, Pyrs only ===
    % useCC = root.info.lyrID == 1 & root.goodind & root.info.uType == 1;
    % plotPhysCompact(root,sess,datT.fpath{i},'goodLyrPyrs',0,root.info.cluster_id(useCC));

    % === Test get_trialSuccessVshuffle ===
    % useF = sessFrst.valTrials(end-49:end);
    % useN = sessLast.valTrials(1:50);
    % lapScoreF(i) = get_trialSuccessVshuffle(sessFrst,datT.rzloc1(i),useF,500);
    % lapScoreN(i) = get_trialSuccessVshuffle(sessLast,datT.rzloc2(i),useN,500);

    % === Export data to binarized spike matrices for aanchal ===
    % handT = import_xldat("D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\clean_final_cohort","data_curation.xlsx");
    % 
    % nShanks = numel(unique(root.info.shankID));
    % 
    % % % === Assign regions for each unit ===
    % roiUnits = zeros(height(root.info),1);
    % for j = 1:nShanks
    %     tmpUnits = root.info.shankID == j-1;
    %     if datT{i,6+j}{1} == 'ca1'
    %         roiUnits(tmpUnits) = 1;
    %     elseif datT{i,6+j}{1} == 'sub'
    %         roiUnits(tmpUnits) = 2;
    %     end
    % end
    % 
    % % Check for hand-removed units
    % if istable(handT)
    %     tInds = handT.Subject == str2num(datT.mouse{i}(end-2:end));
    %     rmCCs  = handT.unit(tInds);
    %     rmInds = false(height(root.info),1);
    %     for j = 1:length(rmCCs)
    %         rmInds(root.info.cluster_id == rmCCs(j)) = true;
    %     end
    % end
    % 
    % % Pick units to use
    % pkExclude = frstHalf.truePk < 1 & lastHalf.truePk < 1;
    % ca1Exclude = roiUnits == 1 & root.info.fr > 10;
    % disp(['CA1 units flagged for FR violation: ' num2str(root.good(ca1Exclude(root.goodind))')])
    % useUnits = root.info.lyrID(root.goodind) == 1 & root.info.fr(root.goodind) > 0.1 ...
    %     & root.info.uType(root.goodind) & ~pkExclude & ~rmInds(root.goodind) & ~ca1Exclude(root.goodind);
    % nCCs = length(root.good);
    % 
    % % % === Send units to Aanchal === 
    % ca1CCs = roiUnits(root.goodind) == 1 & useUnits;
    % subCCs = roiUnits(root.goodind) == 2 & useUnits;
    % 
    % ca1SpkMat = binarizeUnits(root,root.info.cluster_id(ca1CCs),sess);
    % subSpkMat = binarizeUnits(root,root.info.cluster_id(subCCs),sess);
    % 
    % root.regionID = roiUnits;
    % root.ca1CCs   = ca1CCs;
    % root.subCCs   = subCCs;
    % root.rzloc    = [datT.rzloc1(i) datT.rzloc2(i)];
    % 
    % % Save to Z drive
    % saveRoot(root);
    % save([root.name '_spkMat_ca1_pyrs'],'ca1SpkMat')
    % save([root.name '_spkMat_sub_pyrs'],'subSpkMat')

    % === Plot burst metrics and layer units ===
    useCCs = root.good(rootFrst.ca1CCs | rootFrst.subCCs);
    useUnits = rootFrst.ca1CCs | rootFrst.subCCs;
    roiUnits = zeros(size(root.good));
    roiUnits(rootFrst.ca1CCs) = 1;
    roiUnits(rootFrst.subCCs) = 2;
    
    layerUnitF = plot_layerUnits(root,useCCs); legend('off')

    bstFig = plot_burstMetrics(root,sess,useUnits,roiUnits);
    if i == 1; legend('CA1','Sub'); end

    fsave(layerUnitF,[root.name '_useUnitsXlayer'],1,1);
    fsave(bstFig,[root.name '_burstMetrics'],1,1);
    close all

end

% === Test get_trialSuccessVshuffle ===
% output.lapScoreF = lapScoreF;
% output.lapScoreN = lapScoreN;
end