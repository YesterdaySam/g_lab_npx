function [] = combine_rzShiftDat(datT,fname,sdir,sessType)
%% 
% Inputs:
%   datT = a table organizing sessions by mouse and recording day
%   fname = name of the file containing data variables
%   (Deprecated) region = string, e.g. 'ca1' for matching the region ID of each shank in datT
%   sdir = location to save combined data variables
%   sessType = 1 = Fixed RZ; 2 = RZ Shift; 3 = RZ Rand
%
% Outputs:
%   None (variables saved in-function)
% Updated 5/5/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

overwrite = 0;
recID = [];     % [mouseID, recDay, recUnit ID, dist2center, dist2border]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
rgDat = [];     % [recording region for each unit 1 = ca1; 2 = sub; 0 = other
lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
lcMap = [];     % [frstHalf.posfr, lastHalf.posfr];
rpDat = [];     % [frstHalf.ripParticip, frstHalf.ripModBin, lastHalf.ripParticip, lastHalf.ripModBin]
rpRat = [];     % [frstHalf.ripRate, lastHalf.ripRate];
rpMap = [];     % [frstHalf.swrfr, lastHalf.swrfr];
rpMapZ = [];    % [frstHalf.swrz, lastHalf.swrz];
% rpMod = [];     % [frstHalf.ripModBins, lastHalf.ripModBins];
pvStr = [];     % [frstHalf.pvXlap, lastHalf.pvXlap] Struc containing PVC data per lap
frDat = [];     % [frstHalf.standFR, frstHalf.runFR, lastHalf.standFR, lastHalf.runFR]

% recID = [];     % [mouseID, recDay, recUnit ID, dist2center, dist2border]
% useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
% lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
% lcMap = [];     % [frstHalf.posfr, lastHalf.posfr];
% % siStr = [];     % Struc containing binned SI data per 10 laps
% vlDat = [];     % [frstHalf sig., frstHalf slope, frstHalf R, lastHalf sig., lastHalf slope, lastHalf R]
% thDat = [];     % [frstHalf.p, frstHalf.mrl, frstHalf.ang, lastHalf.p, lastHalf.mrl, lastHalf.ang]
% thMap = [];     % [frstHalf.thetafr, lastHalf.thetafr];
% rpDat = [];     % [frstHalf.ripParticip, frstHalf.ripModBin, lastHalf.ripParticip, lastHalf.ripModBin]
% rpRat = [];     % [frstHalf.ripRate, lastHalf.ripRate];
% rpMap = [];     % [frstHalf.swrfr, lastHalf.swrfr];
% rpMapZ = [];    % [frstHalf.swrz, lastHalf.swrz];
% rpMod = [];     % [frstHalf.ripModBins, lastHalf.ripModBins];
% bsDat = [];     % [frstHalf.burstIndex, lastHalf.burstIndex];
% bvDat = [];     % [frstHalf.lckDI, frstHalf.preRZV, lastHalf.lckDI, lastHalf.preRZV
% pvStr = [];     % [frstHalf.pvXlap, lastHalf.pvXlap] Struc containing PVC data per lap

% Load in previously saved data if crashed in midst of run
if overwrite == 0
    try
        cd(sdir)
        physfile = dir('*_phys.mat');
        load(physfile.name)
    catch
    end
else
    ct    = 1;
end

for i = ct:height(datT)

    % === Load data ===
    cd(datT.fpath{i})

    rootfile = dir("*_root.mat");
    load(rootfile.name)
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    epochfile = dir("*_dat.mat");
    load(epochfile.name)
    disp(root.name)

    nShanks = numel(unique(root.info.shankID));

    % % === Assign regions for each unit ===
    roiUnits = zeros(height(root.info),1);
    for j = 1:nShanks
        tmpUnits = root.info.shankID == j-1;
        if datT{i,6+j}{1} == 'ca1'
            roiUnits(tmpUnits) = 1;
        elseif datT{i,6+j}{1} == 'sub'
            roiUnits(tmpUnits) = 2;
        end
    end
    rgDat = [rgDat; roiUnits(root.goodind)];

    pkExclude = frstHalf.truePk < 1 & lastHalf.truePk < 1;
    useUnits = root.info.lyrID(root.goodind) == 1 & root.info.fr(root.goodind) > 0.1 ...
        & root.info.uType(root.goodind) & ~pkExclude;
    nCCs = length(root.good);

    % === Concatenate recording data ===
    useCC = logical([useCC; useUnits]);
    recID = [recID; str2num(datT.mouse{i}(end-2:end))*ones(nCCs,1), ...
        datT.session(i)*ones(nCCs,1), root.good]; %, frstHalf.d2cs, frstHalf.d2bs];

    % % === Concatenate Velocity-FR data ===
    % for j = 1:length(root.good)
    %     vlDat = [vlDat; frstHalf.trueVelMdl(j).p, frstHalf.trueVelMdl(j).b, frstHalf.trueVelMdl(j).r ...
    %         lastHalf.trueVelMdl(j).p, lastHalf.trueVelMdl(j).b, lastHalf.trueVelMdl(j).r];
    % end

    % === Concatenate FR data ===
    frDat = [frDat; frstHalf.frStandRun, lastHalf.frStandRun];

    % % === Concatenate Theta Modulation data ===
    % [thAng1, thMRL1, thP1] = get_thAng(frstHalf.thetastats);
    % [thAng2, thMRL2, thP2] = get_thAng(lastHalf.thetastats);
    % thDat = [thDat; thP1', thMRL1', thAng1', thP2', thMRL2', thAng2'];
    % thMap = [thMap; frstHalf.thetafr lastHalf.thetafr];

    % === Concatenate SI and Peak data ===
    try
        lcDat = [lcDat; frstHalf.sigSI, frstHalf.trueSI, frstHalf.truePk, frstHalf.trueLc, ...
            lastHalf.sigSI, lastHalf.trueSI, lastHalf.truePk, lastHalf.trueLc];
    catch
        lcDat = [lcDat; frstHalf.sig, frstHalf.trueSI, frstHalf.truePk, frstHalf.trueLc, ...
            lastHalf.sig, lastHalf.trueSI, lastHalf.truePk, lastHalf.trueLc];
    end
    lcMap = [lcMap; frstHalf.posfr, lastHalf.posfr];

    % === Concatenate SPWR Modulation data ===
    % for j = 1:length(root.good)
    %     cc = root.good(j);
    %     % frstHalf.ripParticipation(j) = get_RipParticipation(rootFrst,cc,sessFrst,root.ripRef,150);
    %     % lastHalf.ripParticipation(j) = get_RipParticipation(rootLast,cc,sessLast,root.ripRef,150);
    %     % [spwr_frst,~,frstHalf.swrz(j,:)] = plot_frXripple(rootFrst,cc,sessFrst,root.ripRef,150,5,0);
    %     % [spwr_last,~,lastHalf.swrz(j,:)] = plot_frXripple(rootLast,cc,sessLast,root.ripRef,150,5,0);
    %     % frstHalf.ripModBins(j,:) = get_confband(squeeze(frstHalf.shufSPWR(:,j,:)),spwr_frst);
    %     % lastHalf.ripModBins(j,:) = get_confband(squeeze(lastHalf.shufSPWR(:,j,:)),spwr_last);
    %     % frstHalf.ripModBinCt(i) = sum(frstHalf.ripModBins);
    %     % lastHalf.ripModBinCt(i) = sum(lastHalf.ripModBins);
    %     % frstHalf.burstIndex(j) = get_burstIndex(rootFrst,sessFrst,cc);
    %     % lastHalf.burstIndex(j) = get_burstIndex(rootLast,sessLast,cc);
    % end

    rpDat = [rpDat; frstHalf.ripParticipation', lastHalf.ripParticipation'];
    % rpDat = [rpDat; frstHalf.ripParticipation', frstHalf.ripModBinCt', ...
    %     lastHalf.ripParticipation', lastHalf.ripModBinCt'];
    rpRat = [rpRat; frstHalf.ripRate, lastHalf.ripRate];
    % rpMap = [rpMap; frstHalf.swrfr, lastHalf.swrfr];
    rpMapZ = [rpMapZ; frstHalf.swrz lastHalf.swrz];
    % rpMod = [rpMod; frstHalf.ripModBins lastHalf.ripModBins];

    % bsDat = [bsDat; frstHalf.burstIndex', lastHalf.burstIndex'];

    % === Concatenate PV data ===
    try
        bothSIUnits = useUnits & (lastHalf.sigSI <= 0.05 & frstHalf.sigSI <= 0.05);
    catch
        bothSIUnits = useUnits & (lastHalf.sig <= 0.05 & frstHalf.sig <= 0.05);
    end
    bothCA1 = bothSIUnits & roiUnits(root.goodind) == 1;
    bothSub = bothSIUnits & roiUnits(root.goodind) == 2;
    idMat = logical(eye(length(frstHalf.binpos)));

    if sum(bothCA1) > 1 && ~isfield(frstHalf,'pvXlapca1')
        frstHalf.pvOddca1 = (squeeze(mean(frstHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        frstHalf.pvEvnca1 = (squeeze(mean(frstHalf.frMap(2:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        lastHalf.pvOddca1 = (squeeze(mean(lastHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        lastHalf.pvEvnca1 = (squeeze(mean(lastHalf.frMap(2:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        uPVPrePreCA1 = get_pvXtime(frstHalf.posfr,frstHalf.frMap,bothCA1);
        uPVPrePstCA1 = get_pvXtime(lastHalf.posfr,frstHalf.frMap,bothCA1,idMat);
        uPVPstPreCA1 = get_pvXtime(frstHalf.posfr,lastHalf.frMap,bothCA1,idMat);
        uPVPstPstCA1 = get_pvXtime(lastHalf.posfr,lastHalf.frMap,bothCA1);
        frstHalf.pvXlapca1 = [uPVPrePreCA1, uPVPrePstCA1];
        lastHalf.pvXlapca1 = [uPVPstPreCA1, uPVPstPstCA1];
    end

    if sum(bothSub) > 1 && ~isfield(frstHalf,'pvXlapsub')
        frstHalf.pvOddsub = (squeeze(mean(frstHalf.frMap(1:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        frstHalf.pvEvnsub = (squeeze(mean(frstHalf.frMap(2:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        lastHalf.pvOddsub = (squeeze(mean(lastHalf.frMap(1:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        lastHalf.pvEvnsub = (squeeze(mean(lastHalf.frMap(2:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        uPVPrePreSub = get_pvXtime(frstHalf.posfr,frstHalf.frMap,bothSub);
        uPVPrePstSub = get_pvXtime(lastHalf.posfr,frstHalf.frMap,bothSub,idMat);
        uPVPstPreSub = get_pvXtime(frstHalf.posfr,lastHalf.frMap,bothSub,idMat);
        uPVPstPstSub = get_pvXtime(lastHalf.posfr,lastHalf.frMap,bothSub);
        frstHalf.pvXlapsub = [uPVPrePreSub, uPVPrePstSub];
        lastHalf.pvXlapsub = [uPVPstPreSub, uPVPstPstSub];
    end

    if sessType == 2
        if sum(bothCA1) > 1
            pvStr(ct).preBlockPVca1 = frstHalf.pvXlapca1;
            pvStr(ct).pstBlockPVca1 = lastHalf.pvXlapca1;
            pvStr(ct).preEvnca1     = frstHalf.pvEvnca1; % Not actually the pvc, just the frmap used to make it
            pvStr(ct).preOddca1     = frstHalf.pvOddca1;
            pvStr(ct).pstEvnca1     = lastHalf.pvEvnca1;
            pvStr(ct).pstOddca1     = lastHalf.pvOddca1;
        end
        if sum(bothSub) > 1
            pvStr(ct).preBlockPVsub = frstHalf.pvXlapsub;
            pvStr(ct).pstBlockPVsub = lastHalf.pvXlapsub;
            pvStr(ct).preEvnsub     = frstHalf.pvEvnsub; % Not actually the pvc, just the frmap used to make it
            pvStr(ct).preOddsub     = frstHalf.pvOddsub;
            pvStr(ct).pstEvnsub     = lastHalf.pvEvnsub;
            pvStr(ct).pstOddsub     = lastHalf.pvOddsub;
        end
    elseif sessType == 3
        pvStr(ct).preBlockPV = frstHalf.pvXlapRR;
        pvStr(ct).pstBlockPV = lastHalf.pvXlapRR;
        pvStr(ct).preEvn     = frstHalf.pvEvnRR; % Not actually the pvc, just the frmap used to make it
        pvStr(ct).preOdd     = frstHalf.pvOddRR;
        pvStr(ct).pstEvn     = lastHalf.pvEvnRR;
        pvStr(ct).pstOdd     = lastHalf.pvOddRR;
     end

    ct = ct + 1;

    cd(sdir)
    save(fname,'recID','useCC','rgDat','lcDat','lcMap','rpDat','rpRat','rpMap','rpMapZ','pvStr','frDat','ct')
end

cd(sdir)

end