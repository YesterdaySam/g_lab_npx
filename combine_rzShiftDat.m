function [] = combine_rzShiftDat(datT,fname,sdir)

%% 
frDat = [];     % [frstHalf.standFR, frstHalf.runFR, lastHalf.standFR, lastHalf.runFR]
recID = [];     % [mouseID, recDay, recUnit ID, dist2center, dist2border]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
lcMap = [];     % [frstHalf.posfr, lastHalf.posfr];
siStr = [];     % Struc containing binned SI data per 10 laps
vlDat = [];     % [frstHalf sig., frstHalf slope, frstHalf R, lastHalf sig., lastHalf slope, lastHalf R]
thDat = [];     % [frstHalf.p, frstHalf.mrl, frstHalf.ang, lastHalf.p, lastHalf.mrl, lastHalf.ang]
thMap = [];     % [frstHalf.thetafr, lastHalf.thetafr];
rpDat = [];     % [frstHalf.ripParticip, frstHalf.ripModBin, lastHalf.ripParticip, lastHalf.ripModBin]
rpRat = [];     % [frstHalf.ripRate, lastHalf.ripRate];
rpMap = [];     % [frstHalf.swrfr, lastHalf.swrfr];
bvDat = [];     % [frstHalf.lckDI, frstHalf.preRZV, lastHalf.lckDI, lastHalf.preRZV
pvStr = [];     % Struc containing binned SI data per 10 laps

ct    = 1;

for i = 1:height(datT)
    if datT.sess_type(i) ~= 2
        continue
    end

    % === Load data ===
    cd(datT.fpath{i})

    rootfile = dir("*_root.mat");
    load(rootfile.name)
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    epochfile = dir("*_RwdShift_Data2.mat");
    load(epochfile.name)
    disp(root.name)

    nShanks = numel(unique(root.info.shankID));

    % === Identify units to include ===
    roiUnits = zeros(length(root.good),1);
    for j = 1:nShanks
        tmpUnits = root.info.shankID == j-1;
        if datT{i,6+j}{1} == 'ca1'
            roiUnits(tmpUnits) = 0;
        elseif datT{i,6+j}{1} == 'sub'
            roiUnits(tmpUnits) = 1;
        elseif datT{i,6+j}{1} == 'bdr'
            roiUnits(tmpUnits) = 0;
        end
    end

    pkExclude = frstHalf.truePk < 1 & lastHalf.truePk < 1;
    useUnits = root.info.lyrID(root.goodind) == 1 & root.info.fr(root.goodind) > 0.1 ...
        & root.info.uType(root.goodind) & roiUnits(root.goodind) & ~pkExclude;
    nCCs = length(root.good);

    % === Concatenate recording data ===
    useCC = logical([useCC; useUnits]);
    recID = [recID; str2num(datT.mouse{i}(end-2:end))*ones(nCCs,1), ...
        datT.session(i)*ones(nCCs,1), root.good, frstHalf.d2cs, frstHalf.d2bs];

    % === Concatenate Velocity-FR data ===
    for j = 1:length(root.good)
        vlDat = [vlDat; frstHalf.trueVelMdl(j).p, frstHalf.trueVelMdl(j).b, frstHalf.trueVelMdl(j).r ...
            lastHalf.trueVelMdl(j).p, lastHalf.trueVelMdl(j).b, lastHalf.trueVelMdl(j).r];
    end

    % === Concatenate FR data ===
    frDat = [frDat; frstHalf.frStandRun, lastHalf.frStandRun];

    % === Concatenate Theta Modulation data ===
    [thAng1, thMRL1, thP1] = get_thAng(frstHalf.thetastats);
    [thAng2, thMRL2, thP2] = get_thAng(lastHalf.thetastats);
    thDat = [thDat; thP1', thMRL1', thAng1', thP2', thMRL2', thAng2'];
    % for j = 1:length(root.good)
    %     thDat = [thDat; frstHalf.thetastats(j).p, frstHalf.thetastats(j).mrl, frstHalf.thetastats(j).ang, ...
    %         lastHalf.thetastats(j).p, lastHalf.thetastats(j).mrl, lastHalf.thetastats(j).ang];
    % end
    thMap = [thMap; frstHalf.thetafr lastHalf.thetafr];

    % === Concatenate SI and Peak data ===
    lcDat = [lcDat; frstHalf.sig, frstHalf.trueSI, frstHalf.truePk, frstHalf.trueLc, ...
        lastHalf.sig, lastHalf.trueSI, lastHalf.truePk, lastHalf.trueLc];
    lcMap = [lcMap; frstHalf.posfr, lastHalf.posfr];

    % === Concatenate SPWR Modulation data ===
    rpDat = [rpDat; frstHalf.ripParticipation', frstHalf.ripModBinCt', ...
        lastHalf.ripParticipation', lastHalf.ripModBinCt'];
    rpRat = [rpRat; frstHalf.ripRate, lastHalf.ripRate];
    rpMap = [rpMap; frstHalf.swrfr lastHalf.swrfr];

    % === Concatenate SI over 10-trial blocks ===
    siStr(ct).preBlockSI = frstHalf.subEpochSI;
    siStr(ct).pstBlockSI = lastHalf.subEpochSI;

    % === Concatenate PV over laps ===
    pvStr(ct).preBlockPV = frstHalf.pvXlap;
    pvStr(ct).pstBlockPV = lastHalf.pvXlap;
    
    % === Concatenate Behavior data ===
    bvDat(ct).preLckDI = frstHalf.lckDI;
    bvDat(ct).uPreLckDI = frstHalf.ulckDI;
    bvDat(ct).preRZVel = frstHalf.preRZV;
    bvDat(ct).uPreRZVel = frstHalf.uPreRZV;
    bvDat(ct).pstLckDI = lastHalf.lckDI;
    bvDat(ct).uPstLckDI = lastHalf.ulckDI;
    bvDat(ct).pstRZVel = lastHalf.preRZV;
    bvDat(ct).uPstRZVel = lastHalf.uPreRZV;
    
    ct = ct + 1;
end

cd(sdir)

save(fname,'frDat','recID','useCC','siStr','lcDat','lcMap','thDat','thMap','rpDat','rpRat','rpMap','vlDat','bvDat','pvStr')

end