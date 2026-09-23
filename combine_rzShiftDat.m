function [] = combine_rzShiftDat(datT,fname,sdir,sessType,overwrite)
%%
% Inputs:
%   datT = a table organizing sessions by mouse and recording day
%   fname = name of the file containing data variables
%   sdir = location to save combined data variables
%   sessType = 1 = Fixed RZ; 2 = RZ Shift; 3 = RZ Rand
%   overwrite = binary, whether to overwrite or use previously saved vars
%
% Outputs:
%   None (variables saved in-function)
%
% Updated 5/5/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    datT
    fname
    sdir
    sessType
    overwrite = 0;
end

recID = [];     % [mouseID, recDay, recUnit ID, dist2center, dist2border]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
rgDat = [];     % [recording region for each unit 1 = ca1; 2 = sub; 0 = other
lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
lcMap = [];     % [frstHalf.posfr, lastHalf.posfr];
rpDat = [];     % [frstHalf.ripParticip, frstHalf.ripModBin, lastHalf.ripParticip, lastHalf.ripModBin]
rpRat = [];     % [frstHalf.ripRate, frstHalf.uRipDur, frstHalf.pLongSWR, lastHalf.ripRate, lastHalf.uRipDur, lastHalf.pLongSWR];
rpMap = [];     % [frstHalf.swrfr, lastHalf.swrfr];
rpMapZ = [];    % [frstHalf.swrz, lastHalf.swrz];
% rpMod = [];     % [frstHalf.ripModBins, lastHalf.ripModBins];
pvStr = [];     % [frstHalf.pvXlap, lastHalf.pvXlap] Struc containing PVC data per lap
frDat = [];     % [frstHalf.standFR, frstHalf.runFR, lastHalf.standFR, lastHalf.runFR]
bsDat = [];     % [frstHalf.burstIndex, burstISI, burstLen, lastHalf.burstIndex, burstISI, burstLen];
dcDat = [];     % [fxfDecode_errLoc nxnDecode_errLoc nxfDecode_errLoc fxfDecode_errAvg nxnDecode_errAvg nxfDecode_errAvg]
vlDat = [];     % [frstHalf sig., frstHalf slope, frstHalf R, lastHalf sig., lastHalf slope, lastHalf R]

% siStr = [];     % Struc containing binned SI data per 10 laps
% thDat = [];     % [frstHalf.p, frstHalf.mrl, frstHalf.ang, lastHalf.p, lastHalf.mrl, lastHalf.ang]
% thMap = [];     % [frstHalf.thetafr, lastHalf.thetafr];
% rpMod = [];     % [frstHalf.ripModBins, lastHalf.ripModBins];

% Load in previously saved data if crashed in midst of run
if overwrite == 0
    try
        cd(sdir)
        physfile = dir(['*' fname '.mat']);
        load(physfile.name) % saved with ct value
    catch
        ct = 1;
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

    % === Concatenate burst data ===
    bsDat = [bsDat; frstHalf.burstIndex',frstHalf.burstISI',frstHalf.burstLen', ...
        lastHalf.burstIndex',lastHalf.burstISI',lastHalf.burstLen'];
    % % Make a burst metric figure
    % bstFig = plot_burstMetrics(root,sess,useUnits,roiUnits);
    % if ct == 1; legend('CA1','Sub'); end
    % fsave(bstFig,[root.name '_burstMetrics'],1,1);
    % close all

    % % === Concatenate Velocity-FR data ===
    velStats = [];
    for j = 1:length(root.good)
        velStats(j,:) = [frstHalf.trueVelMdl(j).p, frstHalf.trueVelMdl(j).b, frstHalf.trueVelMdl(j).r, ...
            lastHalf.trueVelMdl(j).p, lastHalf.trueVelMdl(j).b, lastHalf.trueVelMdl(j).r];
        vlDat = [vlDat; velStats];
    end

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
    % wlen = 125;
    % if size(frstHalf.shufSPWR,3) == 60
    %     frstHalf = rmfield(frstHalf,{'shufSPWR'});
    %     lastHalf = rmfield(lastHalf,{'shufSPWR'});
    %
    %     frstHalf = get_shufParams(rootFrst,rootFrst.good,sessFrst,frstHalf,250,false,true,false,[],5,wlen);
    %     lastHalf = get_shufParams(rootLast,rootLast.good,sessLast,lastHalf,250,false,true,false,[],5,wlen);
    % end
    % try
    %     frstHalf = rmfield(frstHalf,{'swrz'});
    %     lastHalf = rmfield(lastHalf,{'swrz'});
    % catch
    % end
    % try
    %     frstHalf = rmfield(frstHalf,{'swrfr'});
    %     lastHalf = rmfield(lastHalf,{'swrfr'});
    % catch
    % end
    % try
    %     frstHalf = rmfield(frstHalf,{'ripModBins'});
    %     lastHalf = rmfield(lastHalf,{'ripModBins'});
    % catch
    % end
    % frstHalf = get_unitParams(rootFrst,rootFrst.good,sessFrst,frstHalf,false,false,false,true,false,[],false,0.05,wlen,5);
    % lastHalf = get_unitParams(rootLast,rootLast.good,sessLast,lastHalf,false,false,false,true,false,[],false,0.05,wlen,5);
    %
    % for j = 1:length(root.good)
    %     cc = root.good(j);
    %     frstHalf.ripParticipation(j) = get_RipParticipation(rootFrst,cc,sessFrst,root.ripRef,wlen);
    %     lastHalf.ripParticipation(j) = get_RipParticipation(rootLast,cc,sessLast,root.ripRef,wlen);
    %     % [~,~,frstHalf.swrz(j,:)] = plot_frXripple(rootFrst,cc,sessFrst,root.ripRef,wlen,5,0);
    %     % [~,~,lastHalf.swrz(j,:)] = plot_frXripple(rootLast,cc,sessLast,root.ripRef,wlen,5,0);
    %     frstHalf.ripModBins(j,:) = get_confband(squeeze(frstHalf.shufSPWR(:,j,:)),frstHalf.swrfr(j,:));
    %     lastHalf.ripModBins(j,:) = get_confband(squeeze(lastHalf.shufSPWR(:,j,:)),lastHalf.swrfr(j,:));
    %     frstHalf.ripModBinCt(j) = sum(frstHalf.ripModBins(j,:));
    %     lastHalf.ripModBinCt(j) = sum(lastHalf.ripModBins(j,:));
    % end
    % rlenFrst = (rootFrst.ripStruc(rootFrst.ripRef).ripples(:,3) - rootFrst.ripStruc(rootFrst.ripRef).ripples(:,1)) / rootFrst.fs_lfp *1000;
    % rlenLast = (rootLast.ripStruc(rootLast.ripRef).ripples(:,3) - rootLast.ripStruc(rootLast.ripRef).ripples(:,1)) / rootLast.fs_lfp *1000;
    % frstHalf.pLongSWR = sum(rlenFrst > 100) / length(rlenFrst);
    % lastHalf.pLongSWR = sum(rlenLast > 100) / length(rlenLast);
    % frstHalf.uRipDur = mean(rlenFrst);
    % lastHalf.uRipDur = mean(rlenLast);

    rpDat = [rpDat; frstHalf.ripParticipation', frstHalf.ripModBinCt', ...
        lastHalf.ripParticipation', lastHalf.ripModBinCt'];
    rpRat = [rpRat; frstHalf.ripRate, frstHalf.uRipDur, frstHalf.pLongSWR, ...
        lastHalf.ripRate, lastHalf.uRipDur, lastHalf.pLongSWR];
    % rpMap = [rpMap; frstHalf.swrfr, lastHalf.swrfr];
    rpMapZ = [rpMapZ; frstHalf.swrz lastHalf.swrz];
    % rpMod = [rpMod; frstHalf.ripModBins lastHalf.ripModBins];

    % === Concatenate PV data ===
    try
        bothSIUnits = useUnits & (lastHalf.sigSI <= 0.05 & frstHalf.sigSI <= 0.05);
    catch
        bothSIUnits = useUnits & (lastHalf.sig <= 0.05 & frstHalf.sig <= 0.05);
    end
    bothCA1 = bothSIUnits & roiUnits(root.goodind) == 1;
    bothSub = bothSIUnits & roiUnits(root.goodind) == 2;
    idMat = logical(eye(length(frstHalf.binpos)));
    % idMat = logical(spdiags(ones(1,13),-6:6,37,37));
    % dgMat = logical(spdiags([1 1],[-round(nBins/2) round(nBins/2)],nBins,nBins));
    % dgMat = logical(spdiags(ones(1,13),round(37/2)-6:round(37/2)+6,nBins,nBins) + spdiags(ones(1,13),-round(37/2)-6:-round(37/2)+6,nBins,nBins));

    if sum(bothCA1) > 1 && ~isfield(frstHalf,'pvXlapca1')
        frstHalf.pvOddca1 = (squeeze(mean(frstHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        frstHalf.pvEvnca1 = (squeeze(mean(frstHalf.frMap(2:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        lastHalf.pvOddca1 = (squeeze(mean(lastHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        lastHalf.pvEvnca1 = (squeeze(mean(lastHalf.frMap(2:2:end,:,bothCA1),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothCA1),1,"omitnan")),[],1))';
        uPVPrePreCA1 = get_pvXtime(frstHalf.posfr,frstHalf.frMap,bothCA1,idMat);
        uPVPrePstCA1 = get_pvXtime(lastHalf.posfr,frstHalf.frMap,bothCA1,idMat);
        uPVPstPreCA1 = get_pvXtime(frstHalf.posfr,lastHalf.frMap,bothCA1,idMat);
        uPVPstPstCA1 = get_pvXtime(lastHalf.posfr,lastHalf.frMap,bothCA1,idMat);
        frstHalf.pvXlapca1 = [uPVPrePreCA1, uPVPrePstCA1];
        lastHalf.pvXlapca1 = [uPVPstPreCA1, uPVPstPstCA1];
    end

    if sum(bothSub) > 1 && ~isfield(frstHalf,'pvXlapsub')
        frstHalf.pvOddsub = (squeeze(mean(frstHalf.frMap(1:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        frstHalf.pvEvnsub = (squeeze(mean(frstHalf.frMap(2:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        lastHalf.pvOddsub = (squeeze(mean(lastHalf.frMap(1:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        lastHalf.pvEvnsub = (squeeze(mean(lastHalf.frMap(2:2:end,:,bothSub),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(1:2:end,:,bothSub),1,"omitnan")),[],1))';
        uPVPrePreSub = get_pvXtime(frstHalf.posfr,frstHalf.frMap,bothSub,idMat);
        uPVPrePstSub = get_pvXtime(lastHalf.posfr,frstHalf.frMap,bothSub,idMat);
        uPVPstPreSub = get_pvXtime(frstHalf.posfr,lastHalf.frMap,bothSub,idMat);
        uPVPstPstSub = get_pvXtime(lastHalf.posfr,lastHalf.frMap,bothSub,idMat);
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
            pvStr(ct).preLapMapca1  = frstHalf.frMap(:,:,bothCA1);
            pvStr(ct).pstLapMapca1  = lastHalf.frMap(:,:,bothCA1);
        end
        if sum(bothSub) > 1
            pvStr(ct).preBlockPVsub = frstHalf.pvXlapsub;
            pvStr(ct).pstBlockPVsub = lastHalf.pvXlapsub;
            pvStr(ct).preEvnsub     = frstHalf.pvEvnsub; % Not actually the pvc, just the frmap used to make it
            pvStr(ct).preOddsub     = frstHalf.pvOddsub;
            pvStr(ct).pstEvnsub     = lastHalf.pvEvnsub;
            pvStr(ct).pstOddsub     = lastHalf.pvOddsub;
            pvStr(ct).preLapMapsub  = frstHalf.frMap(:,:,bothSub);
            pvStr(ct).pstLapMapsub  = lastHalf.frMap(:,:,bothSub);
        end
    elseif sessType == 3
        pvStr(ct).preBlockPV = frstHalf.pvXlapRR;
        pvStr(ct).pstBlockPV = lastHalf.pvXlapRR;
        pvStr(ct).preEvn     = frstHalf.pvEvnRR; % Not actually the pvc, just the frmap used to make it
        pvStr(ct).preOdd     = frstHalf.pvOddRR;
        pvStr(ct).pstEvn     = lastHalf.pvEvnRR;
        pvStr(ct).pstOdd     = lastHalf.pvOddRR;
    end

    % Decoding
    useCA1 = useUnits & roiUnits(root.goodind) == 1;
    useSub = useUnits & roiUnits(root.goodind) == 2;
    snameCA1 = [rootFrst.name '_CA1decoders'];
    snameSub = [rootFrst.name '_Subdecoders'];
    % snameSubHi = [rootFrst.name '_Subdecoders_velHi'];
    snameSubLo = [rootFrst.name '_Subdecoders_velLo'];

    splitSlopeSub = prctile(abs(velStats(useSub,2)),50);

    % useSubHi = useSub & abs(velStats(:,2)) >= splitSlopeSub;
    useSubLo = useSub & abs(velStats(:,2)) < splitSlopeSub;

    dcDat(ct).nCA1 = sum(useCA1);
    dcDat(ct).nSub = sum(useSub);
    % dcDat(ct).nSubHi = sum(useSubHi);
    dcDat(ct).nSubLo = sum(useSubLo);

    try
        decodefile = dir("*_Subdecoders_velLo.mat");
        load(decodefile.name)
        [errLocsSub, errMeansSub] = plot_bayesErrHisto(0:0.05:0.95,fxfDecode,nxnDecode,nxfDecode);
    catch
        try
            [errLocsSub, errMeansSub] = decodeRZPlot(useSubLo,snameSubLo,rootFrst,rootLast,sessFrst,sessLast,frstHalf,lastHalf);
        catch
            errLocsSub = nan(1,3); errMeansSub = nan(1,3);
        end
    end

    try
        decodefile = dir("*_CA1decoders.mat");
        load(decodefile.name)
        [errLocsCA1, errMeansCA1] = plot_bayesErrHisto(0:0.05:0.95,fxfDecode,nxnDecode,nxfDecode);
    catch
        try
            [errLocsCA1, errMeansCA1] = decodeRZPlot(useCA1,snameCA1,rootFrst,rootLast,sessFrst,sessLast,frstHalf,lastHalf);
        catch
            errLocsCA1 = nan(1,3); errMeansCA1 = nan(1,3);
        end
    end
    
    try
        decodefile = dir("*_CA1decoders.mat");
        load(decodefile.name)
        [dcDat(ct).ca1_fxf_lapAbsErr,dcDat(ct).ca1_fxf_lapRawErr,dcDat(ct).ca1_fxf_lapErrLoc] = get_decodeErrXlaps(fxfDecode,sessFrst);
        [dcDat(ct).ca1_nxn_lapAbsErr,dcDat(ct).ca1_nxn_lapRawErr,dcDat(ct).ca1_nxn_lapErrLoc] = get_decodeErrXlaps(nxnDecode,sessLast);
        [dcDat(ct).ca1_nxf_lapAbsErr,dcDat(ct).ca1_nxf_lapRawErr,dcDat(ct).ca1_nxf_lapErrLoc] = get_decodeErrXlaps(nxfDecode,sessLast);

        errXlapF = figure; hold on;
        plot(1:sessFrst.nlaps, dcDat(ct).ca1_fxf_lapAbsErr,'k')
        plot((1:sessLast.nlaps) + sessFrst.nlaps, dcDat(ct).ca1_nxn_lapAbsErr,'r')
        plot((1:sessLast.nlaps) + sessFrst.nlaps, dcDat(ct).ca1_nxf_lapAbsErr,'c')
        ylabel('Abs Error (m)'); xlabel('Lap #');
        set(gca,'FontSize',12,'FontName','Arial')
        fsave(errXlapF,[root.name '_ca1decodeErrXLap'],1,0);
    catch
        dcDat(ct).ca1_fxf_lapAbsErr = nan(sessFrst.nlaps,1);
        dcDat(ct).ca1_fxf_lapRawErr = nan(sessFrst.nlaps,1);
        dcDat(ct).ca1_fxf_lapErrLoc = nan(sessFrst.nlaps,1);
        dcDat(ct).ca1_nxn_lapAbsErr = nan(sessLast.nlaps,1);
        dcDat(ct).ca1_nxn_lapRawErr = nan(sessLast.nlaps,1);
        dcDat(ct).ca1_nxn_lapErrLoc = nan(sessLast.nlaps,1);
        dcDat(ct).ca1_nxf_lapAbsErr = nan(sessLast.nlaps,1);
        dcDat(ct).ca1_nxf_lapRawErr = nan(sessLast.nlaps,1);
        dcDat(ct).ca1_nxf_lapErrLoc = nan(sessLast.nlaps,1);
    end

    try
        decodefile = dir("*_Subdecoders_velLo.mat");
        load(decodefile.name)
        [dcDat(ct).sub_fxf_lapAbsErr,dcDat(ct).sub_fxf_lapRawErr,dcDat(ct).sub_fxf_lapErrLoc] = get_decodeErrXlaps(fxfDecode,sessFrst);
        [dcDat(ct).sub_nxn_lapAbsErr,dcDat(ct).sub_nxn_lapRawErr,dcDat(ct).sub_nxn_lapErrLoc] = get_decodeErrXlaps(nxnDecode,sessLast);
        [dcDat(ct).sub_nxf_lapAbsErr,dcDat(ct).sub_nxf_lapRawErr,dcDat(ct).sub_nxf_lapErrLoc] = get_decodeErrXlaps(nxfDecode,sessLast);

        errXlapF = figure; hold on;
        plot(1:sessFrst.nlaps, dcDat(ct).sub_fxf_lapAbsErr,'k')
        plot((1:sessLast.nlaps) + sessFrst.nlaps, dcDat(ct).sub_nxn_lapAbsErr,'r')
        plot((1:sessLast.nlaps) + sessFrst.nlaps, dcDat(ct).sub_nxf_lapAbsErr,'c')
        ylabel('Abs Error (m)'); xlabel('Lap #');
        set(gca,'FontSize',12,'FontName','Arial')
        fsave(errXlapF,[root.name '_subdecodeVelLoErrXLap'],1,0);
    catch
        dcDat(ct).sub_fxf_lapAbsErr = nan(sessFrst.nlaps,1);
        dcDat(ct).sub_fxf_lapRawErr = nan(sessFrst.nlaps,1);
        dcDat(ct).sub_fxf_lapErrLoc = nan(sessFrst.nlaps,1);
        dcDat(ct).sub_nxn_lapAbsErr = nan(sessLast.nlaps,1);
        dcDat(ct).sub_nxn_lapRawErr = nan(sessLast.nlaps,1);
        dcDat(ct).sub_nxn_lapErrLoc = nan(sessLast.nlaps,1);
        dcDat(ct).sub_nxf_lapAbsErr = nan(sessLast.nlaps,1);
        dcDat(ct).sub_nxf_lapRawErr = nan(sessLast.nlaps,1);
        dcDat(ct).sub_nxf_lapErrLoc = nan(sessLast.nlaps,1);
    end

    close all

    dcDat(ct).errLocsCA1 = errLocsCA1;
    dcDat(ct).errMeansCA1 = errMeansCA1;
    dcDat(ct).errLocsSub = errLocsSub;
    dcDat(ct).errMeansSub = errMeansSub;

    ct = ct + 1;

    % save([rootFrst.name '_dat'],'frstHalf','lastHalf','rootFrst','sessFrst','rootLast','sessLast')
    cd(sdir)
    save(fname,'recID','useCC','rgDat','lcDat','lcMap','rpDat','rpRat','rpMap','rpMapZ','pvStr','frDat','bsDat','dcDat','vlDat','ct')
end

cd(sdir)

end