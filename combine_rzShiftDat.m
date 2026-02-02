function [] = combine_rzShiftDat(datT,fname,region,sdir)

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
rpMapZ = [];    % [frstHalf.swrz, lastHalf.swrz];
rpMod = [];     % [frstHalf.ripModBins, lastHalf.ripModBins];
bsDat = [];     % [frstHalf.burstIndex, lastHalf.burstIndex];
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

    % Redundant remake first/last root and sess
    rwdShift = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
    if isfield(sess,'rwdTrials')
        rwdShift = sess.rwdTrials(rwdShift);
    end

    frstHalfInds         = [sess.ind(1) sess.lapend(rwdShift)];
    lastHalfInds         = [sess.lapstt(sess.valTrials(rwdShift+1)) sess.ind(end)];
    [sessFrst, rootFrst] = epochStruc(sess,root,frstHalfInds);
    [sessLast, rootLast] = epochStruc(sess,root,lastHalfInds);

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
    thMap = [thMap; frstHalf.thetafr lastHalf.thetafr];

    % === Concatenate SI and Peak data ===
    lcDat = [lcDat; frstHalf.sig, frstHalf.trueSI, frstHalf.truePk, frstHalf.trueLc, ...
        lastHalf.sig, lastHalf.trueSI, lastHalf.truePk, lastHalf.trueLc];
    lcMap = [lcMap; frstHalf.posfr, lastHalf.posfr];

    % === Concatenate SPWR Modulation data ===
    for j = 1:length(root.good)
        cc = root.good(j);
        % [spwr_frst,~,frstHalf.swrz(j,:)] = plot_frXripple(rootFrst,cc,sessFrst,root.ripRef,150,5,0);
        % [spwr_last,~,lastHalf.swrz(j,:)] = plot_frXripple(rootLast,cc,sessLast,root.ripRef,150,5,0);
        % frstHalf.ripModBins(j,:) = get_confband(squeeze(frstHalf.shufSPWR(:,j,:)),spwr_frst);
        % lastHalf.ripModBins(j,:) = get_confband(squeeze(lastHalf.shufSPWR(:,j,:)),spwr_last);
        frstHalf.burstIndex(j) = get_burstIndex(rootFrst,sessFrst,cc);
        lastHalf.burstIndex(j) = get_burstIndex(rootLast,sessLast,cc);
    end

    rpDat = [rpDat; frstHalf.ripParticipation', frstHalf.ripModBinCt', ...
        lastHalf.ripParticipation', lastHalf.ripModBinCt'];
    rpRat = [rpRat; frstHalf.ripRate, lastHalf.ripRate];
    rpMap = [rpMap; frstHalf.swrfr lastHalf.swrfr];
    rpMapZ = [rpMapZ; frstHalf.swrz lastHalf.swrz];
    rpMod = [rpMod; frstHalf.ripModBins lastHalf.ripModBins];

    bsDat = [bsDat; frstHalf.burstIndex', lastHalf.burstIndex'];

    % === Concatenate SI over 10-trial blocks ===
    siStr(ct).preBlockSI = frstHalf.subEpochSI;
    siStr(ct).pstBlockSI = lastHalf.subEpochSI;

    % === Concatenate PV over laps ===
    pvStr(ct).preBlockPV = frstHalf.pvXlap;
    pvStr(ct).pstBlockPV = lastHalf.pvXlap;
    
    % === Concatenate Behavior data ===
    [~,~,preLckMap] = plot_lickpos(sessFrst,0.03,0);
    [~,~,pstLckMap] = plot_lickpos(sessLast,0.03,0);

    bvDat(ct).preLckDI    = frstHalf.lckDI;
    bvDat(ct).uPreLckDI   = frstHalf.ulckDI;
    bvDat(ct).preRZVel    = frstHalf.preRZV;
    bvDat(ct).uPreRZVel   = frstHalf.uPreRZV;
    bvDat(ct).pstLckDI    = lastHalf.lckDI;
    bvDat(ct).uPstLckDI   = lastHalf.ulckDI;
    bvDat(ct).pstRZVel    = lastHalf.preRZV;
    bvDat(ct).uPstRZVel   = lastHalf.uPreRZV;
    bvDat(ct).uPreRZVel20 = frstHalf.uPreRZV20;
    bvDat(ct).uPreLckDI20 = frstHalf.ulckDI20;
    bvDat(ct).uPstRZVel20 = lastHalf.uPreRZV20;
    bvDat(ct).uPstLckDI20 = lastHalf.ulckDI20;
    bvDat(ct).preVMap     = frstHalf.spVelMap;
    bvDat(ct).pstVMap     = lastHalf.spVelMap;
    bvDat(ct).uPreVF20    = mean(frstHalf.spVelMap(1:20,:),'omitnan');
    bvDat(ct).uPreVL20    = mean(frstHalf.spVelMap(end-19:end,:),'omitnan');
    bvDat(ct).uPstVF20    = mean(lastHalf.spVelMap(1:20,:),'omitnan');
    bvDat(ct).uPstVL20    = mean(lastHalf.spVelMap(end-19:end,:),'omitnan');
    bvDat(ct).preLMap     = preLckMap;
    bvDat(ct).pstLMap     = pstLckMap;
    bvDat(ct).uPreLF20    = mean(preLckMap(1:20,:),'omitnan');
    bvDat(ct).uPreLL20    = mean(preLckMap(end-19:end,:),'omitnan');
    bvDat(ct).uPstLF20    = mean(pstLckMap(1:20,:),'omitnan');
    bvDat(ct).uPstLL20    = mean(pstLckMap(end-19:end,:),'omitnan');

    frstHalf.vCorXLap = corr(frstHalf.spVelMap', mean(frstHalf.spVelMap(end-29:end,:))','rows','complete');
    try
        lastHalf.vCorXLap = corr(lastHalf.spVelMap', mean(lastHalf.spVelMap(31:60,:))','rows','complete');
    catch
        lastHalf.vCorXLap = corr(lastHalf.spVelMap', mean(lastHalf.spVelMap(end-29:end,:))','rows','complete');
    end

    bvDat(ct).vCorPre     = frstHalf.vCorXLap;
    bvDat(ct).vCorPst     = lastHalf.vCorXLap;
    
    ct = ct + 1;

    save([root.name '_RwdShift_Data2'],'frstHalf','lastHalf');
end

cd(sdir)

save(fname,'frDat','recID','useCC','siStr','lcDat','lcMap','thDat','thMap','rpDat','rpRat','rpMap','rpMapZ','rpMod','bsDat','vlDat','bvDat','pvStr')

end