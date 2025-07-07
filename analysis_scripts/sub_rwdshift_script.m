
spath = 'D:\Data\Kelton\analyses\KW049\KW049_06172025_rec_D2_RLat2';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)

rwdShift = sess.valTrials(find(diff(sess.pos(sess.rwdind)) > 0.1,1)+1);   % Find lap of reward shift
nUnits = length(root.good);
saveFlag = 1;

dbnsz = 0.05;
histoBnsz = 5;
wlen = 150;
ripRef = root.ripRef;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm

%% Epoch into pre/post reward shift for reward shift sessions

if ~isempty(rwdShift)
    frstHalfInds         = [sess.ind(1) sess.lapend(rwdShift-1)];
    lastHalfInds         = [sess.lapstt(rwdShift) sess.ind(end)];
    [sessFrst, rootFrst] = epochStruc(sess,root,frstHalfInds);
    [sessLast, rootLast] = epochStruc(sess,root,lastHalfInds);
end

%% Plot example pre/post
cc = 95;

rzPosFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,1);
rzVelFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,2);
rzThMFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,3);
rzRwdFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,4);

if saveFlag
    tblind = find(root.info.cluster_id == cc);
    sbase = ['unit_' num2str(cc) '_shank_' num2str(root.info.shankID(tblind)) '_rwdshift_'];
    saveas(rzPosFig,[sbase 'spatial'], 'png')
    saveas(rzVelFig,[sbase 'velocity'],'png')
    saveas(rzThMFig,[sbase 'thetaMod'],'png')
    saveas(rzRwdFig,[sbase 'rewardT'], 'png')
end

%% Get real unit and epoch parameters

clear frstHalf lastHalf
frstHalf.name   = root.name;
lastHalf.name   = root.name;
frstHalf.trueSI = zeros(size(root.good));
frstHalf.truePk = zeros(size(root.good));
frstHalf.trueLc = zeros(size(root.good));
lastHalf.trueSI = zeros(size(root.good));
lastHalf.truePk = zeros(size(root.good));
lastHalf.trueLc = zeros(size(root.good));

for i = 1:nUnits
    cc = root.good(i);
    lfpInd = root.info.shankID(root.info.cluster_id == cc)+1;   % Account for 0-indexing
    [frstHalf.trueSI(i),~,frstHalf.truePk(i),frstHalf.trueLc(i),~,~,~,frstHalf.binedges] = get_SI(rootFrst,cc,sessFrst,dbnsz);
    [lastHalf.trueSI(i),~,lastHalf.truePk(i),lastHalf.trueLc(i),~,~,~,lastHalf.binedges] = get_SI(rootLast,cc,sessLast,dbnsz);
    frstHalf.frStandRun(i,:) = get_frStandVRun(rootFrst,cc,sessFrst);
    lastHalf.frStandRun(i,:) = get_frStandVRun(rootLast,cc,sessLast);
    [~,~,frstHalf.trueVelMdl(i)] = plot_frXvel(rootFrst,cc,sessFrst,2,0);
    [~,~,lastHalf.trueVelMdl(i)] = plot_frXvel(rootLast,cc,sessLast,2,0);
    frstHalf.thetastats(i) = plot_thetaMod(rootFrst,cc,lfpInd,2*pi/36,0);
    lastHalf.thetastats(i) = plot_thetaMod(rootLast,cc,lfpInd,2*pi/36,0);
end

frstHalf = subEpoch(sessFrst,rootFrst,frstHalf);
lastHalf = subEpoch(sessLast,rootLast,lastHalf);

frstHalf.ripRate = size(rootFrst.ripStruc(ripRef).ripples,1) / (sum(not(sessFrst.runInds)) / sess.samprate);    %Normalize based on standing periods
lastHalf.ripRate = size(rootLast.ripStruc(ripRef).ripples,1) / (sum(not(sessLast.runInds)) / sess.samprate);
frstHalf.binpos = frstHalf.binedges(1:end-1)+0.5*dbnsz;
lastHalf.binpos = lastHalf.binedges(1:end-1)+0.5*dbnsz;

%% Find place cells and SPW-R modulation pre and post shift
% Methods: 
% Kitanishi et al. 2021: Exceed 99th percentile of SI shuffle
% Grienberger & Magee 2022:
%   Field: contiguous area with >= 20% peak FR
%   Induction: First lap with >3SD of in-field FR, above out-field noise
%              and activity in 2/5 subsequent laps also >3SD
%   Stability: Significant >3SD activity in 30% of post-induction laps

tic
nShufs = 250;

frstHalf.shufSI = zeros(nUnits,nShufs);
lastHalf.shufSI = zeros(nUnits,nShufs);
frstHalf.shufSPWR = zeros(nShufs,nUnits,length(-wlen:histoBnsz:wlen)-1);
lastHalf.shufSPWR = zeros(nShufs,nUnits,length(-wlen:histoBnsz:wlen)-1);

for j = 1:nShufs
    [shiftFrst,shiftSessFrst,shiftIndFrst] = shiftTrain(rootFrst,sessFrst);
    [shiftLast,shiftSessLast,shiftIndLast] = shiftTrain(rootLast,sessLast);

    if mod(j,50) == 0
        disp(['Shuffle # ' num2str(j)])
        toc 
    end

    for i = 1:nUnits
        cc = root.good(i);

        [frstHalf.shufSI(i,j)] = get_SI(shiftFrst,cc,shiftSessFrst);
        [lastHalf.shufSI(i,j)] = get_SI(shiftLast,cc,shiftSessLast);

        frstHalf.shufSPWR(j,i,:) = plot_frXripple(shiftFrst,cc,shiftSessFrst,ripRef,wlen,histoBnsz,0);
        lastHalf.shufSPWR(j,i,:) = plot_frXripple(shiftLast,cc,shiftSessLast,ripRef,wlen,histoBnsz,0);
    end
end

frstHalf.sig = sum(frstHalf.shufSI > frstHalf.trueSI,2) / nShufs;
lastHalf.sig = sum(lastHalf.shufSI > lastHalf.trueSI,2) / nShufs;

disp(['First half SI p <= 0.05 ' num2str(sum(frstHalf.sig <= 0.05)) ' of ' num2str(nUnits) ' units'])
disp(['Last half SI p <= 0.05 ' num2str(sum(lastHalf.sig <= 0.05)) ' of ' num2str(nUnits) ' units'])
disp(['First and last half SI p <= 0.05 ' num2str(sum(frstHalf.sig <= 0.05 & lastHalf.sig <= 0.05)) ' of ' num2str(nUnits) ' units'])

%% SPWRs pre and post

for i = 1:nUnits
    cc = root.good(i);
    spwr_frst = plot_frXripple(rootFrst,cc,sessFrst,ripRef,wlen,histoBnsz,0);
    frstHalf.ripParticipation(i) = get_RipParticipation(rootFrst,cc,sessFrst,ripRef,wlen);
    spwr_last = plot_frXripple(rootLast,cc,sessLast,ripRef,wlen,histoBnsz,0);
    lastHalf.ripParticipation(i) = get_RipParticipation(rootLast,cc,sessLast,ripRef,wlen);
    p_frst = get_confband(squeeze(frstHalf.shufSPWR(:,i,:)),spwr_frst);
    p_last = get_confband(squeeze(lastHalf.shufSPWR(:,i,:)),spwr_last);
    frstHalf.ripModBinCt(i) = sum(p_frst);
    lastHalf.ripModBinCt(i) = sum(p_last);
end

disp(['First half SPWR Modulated p <= 0.05 ' num2str(sum(frstHalf.ripModBinCt > 1)) ' of ' num2str(nUnits) ' units'])
disp(['Last half SPWR Modulated p <= 0.05 ' num2str(sum(lastHalf.ripModBinCt > 1)) ' of ' num2str(nUnits) ' units'])
disp(['First and last half SPWR Modulated p <= 0.05 ' num2str(sum(frstHalf.ripModBinCt > 1 & lastHalf.ripModBinCt > 1)) ' of ' num2str(nUnits) ' units'])

%% Save for later

if saveFlag
    save([root.name '_RwdShift_Data2'],'frstHalf','lastHalf')
end

%% SPWRs over time

figure; hold on
spwrs_trial = histcounts(root.ripStruc(ripRef).ripples(:,2),[sess.lapstt(sess.valTrials(1):sess.valTrials(end)); sess.lapend(sess.valTrials(end))]);
plot(spwrs_trial,'k')
plot([rwdShift rwdShift],[0 max(spwrs_trial)],'k--')
xlabel('Trial #'); ylabel('SPW-R count')
set(gca,'FontSize',12,'FontName','Arial')

figure; hold on
binedges = sess.ts(1):60:sess.ts(end);
spwrs_time = histcounts(sess.ts(root.ripStruc(ripRef).ripples(:,2)),binedges);
plot(binedges(1:end-1),spwrs_time,'k')
plot([sess.ts(sess.lapstt(rwdShift)) sess.ts(sess.lapstt(rwdShift))],[0 max(spwrs_time)],'k--')
xlabel('Time (s)'); ylabel('SPW-R count')
set(gca,'FontSize',12,'FontName','Arial')

%% Descriptive statistics and graphs
% root.good(sigFrst <= 0.05 & sigLast <= 0.05 & root.info.fr(root.goodind)' > 0.1 & root.info.lyrID(root.goodind)' == 1)

lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = lyrUnits & hiFRUnits & root.info.uType(root.goodind);
siUnits = useUnits & (frstHalf.sig <= 0.05 | lastHalf.sig <= 0.05);
frstSIUnits = useUnits & (frstHalf.sig <= 0.05);
lastSIUnits = useUnits & (lastHalf.sig <= 0.05);
bothSIUnits = useUnits & (lastHalf.sig <= 0.05 & frstHalf.sig <= 0.05);

nUseUnits = sum(useUnits);

xcoords = ones(nUseUnits,1);
shiftbins = [-fliplr(frstHalf.binedges) frstHalf.binedges(2:end)];

% Firing rate stand vs run
[~,ps.standFR_firstlast] = ttest(frstHalf.frStandRun(useUnits,1),lastHalf.frStandRun(useUnits,1));
[~,ps.runFR_firstlast] = ttest(frstHalf.frStandRun(useUnits,2),lastHalf.frStandRun(useUnits,2));

bardat = [mean(frstHalf.frStandRun(useUnits,1)), mean(frstHalf.frStandRun(useUnits,2)); mean(lastHalf.frStandRun(useUnits,1)), mean(lastHalf.frStandRun(useUnits,2))];

FRStandRunFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
b1 = bar(bardat');
b1(1).FaceColor = [0.5 0.5 1];
b1(2).FaceColor = [0.75 0.75 1];
plot(xcoords-0.15,frstHalf.frStandRun(useUnits,1),'ko')
plot(xcoords+0.15,lastHalf.frStandRun(useUnits,1),'ko')
plot(xcoords+1-0.15,frstHalf.frStandRun(useUnits,2),'ko')
plot(xcoords+1+0.15,lastHalf.frStandRun(useUnits,2),'ko')
% plot([frstHalf.frStandRun(useUnits,1) frstHalf.frStandRun(useUnits,2)]','k-o')
xlim([0.5 2.5])
xticks(1:2); xticklabels({'Standing', 'Running'})
ylabel('Firing Rate (Hz)')
legend({'Familiar','Novel'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

% Spatial Information
[~,ps.SI_firstlast] = ttest(frstHalf.trueSI(siUnits),lastHalf.trueSI(siUnits));

bardat = [mean(frstHalf.trueSI(siUnits)); mean(lastHalf.trueSI(siUnits))];

lcBothSIFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
b2 = bar(bardat,'FaceColor','flat','BarWidth',0.5);
b2.CData(1,:) = [0.5 0.5 1];
b2.CData(2,:) = [0.75 0.75 1];
plot([frstHalf.trueSI(siUnits); lastHalf.trueSI(siUnits)],'k-o')
xlim([0.5 2.5]); ylim([-0.01 inf])
ylabel('Spatial Information (Bits/spike)')
xticks(1:2); xticklabels({'Familiar', 'Novel'})
set(gca,'FontSize',12,'FontName','Arial')

% % Peak location delta (positive indicates field shifted forward in last half)
% deltas_Lc = lastHalf.binpos(lastHalf.trueLc) - frstHalf.binpos(frstHalf.trueLc);
% delta_counts = histcounts(deltas_Lc(siUnits),shiftbins);
% 
% peakDeltaFig = figure; hold on
% bar((shiftbins(1:end-1)+0.5*dbnsz)*100, delta_counts./sum(delta_counts),'FaceColor',[0.55, 0.45, 1])
% ylabel('Probability')
% xlabel('Peak FR Shift Novel - Familiar (cm)')
% xlim([shiftbins(1) shiftbins(end)]*100)
% set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    sbase = root.name;
    saveas(FRStandRunFig,[sbase '_RwdShift_FRStandRun.png'])
    saveas(lcBothSIFig,[sbase '_RwdShift_SI.png'])
    % saveas(peakDeltaFig,[sbase '_RwdShift_PeakShiftRZ.png'])
end

%% Peak Location distribution
fieldDistroFrst = histcounts(frstHalf.binpos(frstHalf.trueLc(siUnits)),frstHalf.binedges);
fieldDistroLast = histcounts(lastHalf.binpos(lastHalf.trueLc(siUnits)),lastHalf.binedges);

peakDistroFig = figure; hold on
plot(frstHalf.binpos*100,fieldDistroFrst./sum(fieldDistroFrst),'Color',[0.5 0.5 1],'LineWidth',2);
plot(lastHalf.binpos*100,fieldDistroLast./sum(fieldDistroLast),'Color',[0.75 0.75 1],'LineWidth',2);
plot([r1pos r1pos]*100,[0 0.15],'--','Color',[0.5 0.5 1])
plot([r2pos r2pos]*100,[0 0.15],'--','Color',[0.75 0.75 1])
xlabel('Track Position (cm)')
ylabel('Probability of field peak')
legend({'Familiar','Novel','RZ1','RZ2'},'Location','northeast')
set(gca,'FontSize',12,'FontName','Arial')

% Align field distribution around reward
trackLen = max(frstHalf.binedges);
shiftR1 = trackLen/2 - r1pos;
fieldAlignR1 = mod(frstHalf.binpos(frstHalf.trueLc(siUnits))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
shiftR2 = trackLen/2 - r2pos;
fieldAlignR2 = mod(lastHalf.binpos(lastHalf.trueLc(siUnits))+shiftR2,trackLen)+0.5*dbnsz;

fieldDistroFrst_Align = histcounts(fieldAlignR1,frstHalf.binedges);
fieldDistroLast_Align = histcounts(fieldAlignR2,lastHalf.binedges);

rzDistroFig = figure; hold on
plot((frstHalf.binpos-trackLen/2)*100,fieldDistroFrst_Align./sum(fieldDistroFrst_Align),'Color',[0.5 0.5 1],'LineWidth',2);
plot((lastHalf.binpos-trackLen/2)*100,fieldDistroLast_Align./sum(fieldDistroLast_Align),'Color',[0.75 0.75 1],'LineWidth',2);
plot([0 0]*100,[0 0.15],'k--')
ylabel('Probability of field peak')
xlabel('Distance to RZ (cm)')
legend({'Familiar','Novel','RZ'},'Location','northeast')
set(gca,'FontSize',12,'FontName','Arial')

% Compare shift in PF peak relative to reward per unit
deltaField_RZAlign = fieldAlignR2 - fieldAlignR1;
circAlignNeg = deltaField_RZAlign < -trackLen/2;
circAlignPos = deltaField_RZAlign > trackLen/2;
deltaField_RZAlign(circAlignNeg) = deltaField_RZAlign(circAlignNeg) + trackLen;     % When new field forward-shifts
deltaField_RZAlign(circAlignPos) = deltaField_RZAlign(circAlignPos) - trackLen;     % When new field back-shifts

deltaField_Distro = histcounts(deltaField_RZAlign,shiftbins);

fieldShiftRZFig = figure; hold on
bar((shiftbins(1:end-1)+0.5*dbnsz)*100,deltaField_Distro/sum(deltaField_Distro),'FaceColor',[0.25 0.15 1])
plot([0 0]*100,[0 0.15],'k--')
xlabel('\Delta RZ-aligned Novel - Familiar (cm)')
ylabel('Probability')
xlim([-trackLen/2-dbnsz trackLen/2+dbnsz]*100)
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    sbase = root.name;
    saveas(peakDistroFig,[sbase '_RwdShift_PeakDistro.png'])
    saveas(rzDistroFig,[sbase '_RwdShift_PeakDistroRZ.png'])
    saveas(fieldShiftRZFig,[sbase '_RwdShift_PeakShiftRZ.png'])
end

%% Waterfall plots
frstHalf.rwdBin = find(frstHalf.binpos > r1pos,1);
lastHalf.rwdBin = find(lastHalf.binpos > r2pos,1);

[frstRZ_SI_Wfl_fig] = plot_unitsXpos(rootFrst,sessFrst,root.good(frstSIUnits));
plot([frstHalf.rwdBin frstHalf.rwdBin],[0 length(frstSIUnits)+1],'r--','LineWidth',2)
title('Familiar RZ')
lastRZ_SI_Wfl_fig = plot_unitsXpos(rootLast,sessLast,root.good(lastSIUnits));
plot([lastHalf.rwdBin lastHalf.rwdBin],[0 length(lastSIUnits)+1],'r--','LineWidth',2)
title('Novel RZ')

[frstRZ_bothSI_Wfl_fig,~,tmpsort] = plot_unitsXpos(rootFrst,sessFrst,root.good(bothSIUnits));
plot([frstHalf.rwdBin frstHalf.rwdBin],[0 length(bothSIUnits)+1],'r--','LineWidth',2)
title('Familiar RZ - Sorted by Familiar')
lastRZ_bothSI_Wfl_fig = plot_unitsXpos(rootLast,sessLast,root.good(bothSIUnits));
plot([lastHalf.rwdBin lastHalf.rwdBin],[0 length(bothSIUnits)+1],'r--','LineWidth',2)
title('Novel RZ - Sorted by Novel')
lastRZ_bothSI_Wfl_sort_fig = plot_unitsXpos(rootLast,sessLast,root.good(bothSIUnits),tmpsort);
plot([lastHalf.rwdBin lastHalf.rwdBin],[0 length(bothSIUnits)+1],'r--','LineWidth',2)
title('Novel RZ - Sorted by Familiar')

if saveFlag
    sbase = root.name;
    saveas(lastRZ_SI_Wfl_fig,[sbase '_RwdShift_waterfall_frstSI.png'])
    saveas(frstRZ_SI_Wfl_fig,[sbase '_RwdShift_waterfall_lastSI.png'])
    saveas(frstRZ_bothSI_Wfl_fig,[sbase '_RwdShift_waterfall_bothSI_frst.png'])
    saveas(lastRZ_bothSI_Wfl_fig,[sbase '_RwdShift_waterfall_bothSI_last.png'])
    saveas(lastRZ_bothSI_Wfl_sort_fig,[sbase '_RwdShift_waterfall_bothSI_last_sortfrst.png'])
end

%% Combine Sessions/Animals

parentDir = "D:\Data\Kelton\analyses\group_analyses"; 
subDatT = import_xldat(parentDir,"dat_include.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift';
cd(groupSDir) 

cleanInds = subDatT.include == 0;
subDatT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;

dbnsz = 0.05;
histoBnsz = 5;
wlen = 150;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

%% 
frDat = [];     % [frstHalf.standFR, frstHalf.runFR, lastHalf.standFR, lastHalf.runFR]
recID = [];     % [mouseID, recDay, recUnit ID]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
siStr = [];     % Struc containing binned SI data per 10 laps
vlDat = [];     % [frstHalf sig., frstHalf slope, frstHalf R, lastHalf sig., lastHalf slope, lastHalf R]
thDat = [];     % [frstHalf.p, frstHalf.mrl, frstHalf.ang, lastHalf.p, lastHalf.mrl, lastHalf.ang]
rpDat = [];     % [frstHalf.ripParticip, frstHalf.ripModBin, lastHalf.ripParticip, lastHalf.ripModBin]
rpRat = [];     % [frstHalf.ripRate, lastHalf.ripRate];

clear ps stats

ct    = 1;

for i = 1:height(subDatT)
    if subDatT.sess_type(i) ~= 2
        continue
    end

    % === Load data ===
    cd(subDatT.fpath{i})

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
        if subDatT{i,6+j}{1} == 'ca1'
            roiUnits(tmpUnits) = 0;
        elseif subDatT{i,6+j}{1} == 'sub'
            roiUnits(tmpUnits) = 1;
        end
    end

    useUnits = root.info.lyrID(root.goodind) == 1 & root.info.fr(root.goodind) > 0.1 & root.info.uType(root.goodind) & roiUnits(root.goodind);
    nCCs = length(root.good);

    % === Concatenate recording data ===
    useCC = logical([useCC; useUnits]);
    recID = [recID; str2num(subDatT.mouse{i}(end-2:end))*ones(nCCs,1), subDatT.session(i)*ones(nCCs,1), root.good];

    % === Concatenate Velocity data ===
    for j = 1:length(root.good)
        vlDat = [vlDat; frstHalf.trueVelMdl(j).p, frstHalf.trueVelMdl(j).b, frstHalf.trueVelMdl(j).r...
            lastHalf.trueVelMdl(j).p, lastHalf.trueVelMdl(j).b, lastHalf.trueVelMdl(j).r];
    end

    % === Concatenate FR data ===
    frDat = [frDat; frstHalf.frStandRun, lastHalf.frStandRun];

    % === Concatenate Theta Modulation data ===
    for j = 1:length(root.good)
        thDat = [thDat; frstHalf.thetastats(j).p, frstHalf.thetastats(j).mrl, frstHalf.thetastats(j).ang,...
            lastHalf.thetastats(j).p, lastHalf.thetastats(j).mrl, lastHalf.thetastats(j).ang];
    end

    % === Concatenate SI and Peak data ===
    lcDat = [lcDat; frstHalf.sig, frstHalf.trueSI, frstHalf.truePk, frstHalf.trueLc,...
        lastHalf.sig, lastHalf.trueSI, lastHalf.truePk, lastHalf.trueLc];

    % === Concatenate SPWR Modulation data ===
    rpDat = [rpDat; frstHalf.ripParticipation', frstHalf.ripModBinCt',...
        lastHalf.ripParticipation', lastHalf.ripModBinCt'];
    rpRat = [rpRat; frstHalf.ripRate, lastHalf.ripRate];

    % === Concatenate SI over 10-trial blocks ===
    siStr(ct).preBlockSI = frstHalf.subEpochSI;
    siStr(ct).pstBlockSI = lastHalf.subEpochSI;

    ct = ct + 1;
end

cd(groupSDir)

save('subRwdShift_data2','frDat','recID','useCC','siStr','lcDat','thDat','rpDat','rpRat','vlDat')

%% Velocity coding
% vlDat: 1&4 = sig.; 2&5 = Slope/B; 3&6 = Corr/R
vlFrstID = useCC & vlDat(:,1) <= 0.05;
vlLastID = useCC & vlDat(:,4) <= 0.05;
vlBothID = vlFrstID & vlLastID;

velPieFig = prepostPie(vlFrstID,vlLastID,useCC);
title("Velocity-modulated units");

vlDltB  = vlDat(vlBothID,5) - vlDat(vlBothID,2);
vlDltR  = vlDat(vlBothID,6) - vlDat(vlBothID,3);

[~,ps.vl_DtB_both,~,stats.vl_DtB_both] = ttest(vlDat(vlBothID,2), vlDat(vlBothID,5));
[~,ps.vl_PPB_eith,~,stats.vl_PPB_eith] = ttest2(vlDat(vlFrstID,2),vlDat(vlLastID,5));
% [~,ps.vl_PPB_only,~,stats.vl_PPB_only] = ttest2(vlDat(vlFrstID,2),vlDat(vlLastID,5));
[~,ps.vl_DtR_both,~,stats.vl_DtR_both] = ttest(vlDat(vlBothID,3), vlDat(vlBothID,6));
[~,ps.vl_PPR_eith,~,stats.vl_PPR_eith] = ttest2(vlDat(vlFrstID,3),vlDat(vlLastID,6));

binedges = -1:0.05:1;
velDltBFig = plotDeltaHisto(vlDat(vlBothID,2),vlDat(vlBothID,5),binedges);
xlabel('\Delta Slope (Novel - Familiar)')
title("Significant units pre & post");

binedges = 0:0.02:1;
velDltRFig = plotDeltaHisto(vlDat(vlBothID,3),vlDat(vlBothID,6),binedges);
ylabel('Probability'); xlabel('\Delta R (Novel - Familiar)')
title("Significant units pre & post");

if saveFlag
    saveas(velPieFig,'subRwdShift_vel_pie','png')
    saveas(velDltBFig,'subRwdShift_vel_deltaB','png')
    saveas(velDltRFig,'subRwdShift_vel_deltaR','png')
end

%% Firing Rate stand vs run
% frDat: 1&3 = standing; 2&4 = running

nUseCC = sum(useCC);

[~,ps.fr_PPStnd_all,~,stats.fr_PPStnd_all] = ttest(frDat(useCC,1),frDat(useCC,3));
[~,ps.fr_PPRunn_all,~,stats.fr_PPRunn_all] = ttest(frDat(useCC,2),frDat(useCC,4));

bardat = [mean(frDat(useCC,1)), mean(frDat(useCC,2)); mean(frDat(useCC,3)), mean(frDat(useCC,4))];
xcoords = ones(sum(useCC),1);

FRStandRunFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
b1 = bar(bardat','FaceColor','flat');
b1(1).CData = vColors2(1,:);
b1(2).CData = vColors2(2,:);
errorbar([0.85 1.15 1.85 2.15],reshape(bardat,4,1),std(frDat)./sqrt(nUseCC),'k.')
% plot(xcoords-0.15,frDat(useCC,1),'ko')
% plot(xcoords+0.15,frDat(useCC,3),'ko')
% plot(xcoords+1-0.15,frDat(useCC,2),'ko')
% plot(xcoords+1+0.15,frDat(useCC,4),'ko')
xlim([0.5 2.5])
xticks(1:2); xticklabels({'Standing', 'Running'})
ylabel('Firing Rate (Hz)')
legend({'Familiar','Novel'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(FRStandRunFig,'subRwdShift_fr_bar','png')
end

%% Theta Modulation
% thDat: 1&4 = sig.; 2&5 = MRL; 3&6 = Angle
thFrstID = useCC & thDat(:,1) <= 0.05;
thLastID = useCC & thDat(:,4) <= 0.05;
thBothID = thFrstID & thLastID;
nThBoth = sum(thBothID);

thetaPie = prepostPie(thFrstID,thLastID,useCC);
title("Sig. Theta-modulated units");

[~,ps.th_PPM_both,~,stats.th_PPM_both] = ttest(thDat(thBothID,2),thDat(thBothID,5));
[~,ps.th_PPA_both,~,stats.th_PPA_both] = ttest(thDat(thBothID,3),thDat(thBothID,6));

thMRLFig = plotBar2(thDat(thBothID,2),thDat(thBothID,5));
ylabel('Theta Mean Resultant Length')

binedges = -0.1:0.01:0.1;
thDltMFig = plotDeltaHisto(thDat(thBothID,2),thDat(thBothID,5),binedges);
xlabel('\Delta Theta MRL (Novel - Familiar)')
title("Significant units pre & post");

% 0/360 = trough, 180/540 = peak
thAngFig = plotBar2(rad2deg(thDat(thBothID,3)),rad2deg(thDat(thBothID,6)));
ylabel('Theta Angle (degree)'); ylim([-180 180])

binedges = rad2deg(-pi/2:pi/36:pi/2);
thDltAFig = plotDeltaHisto(rad2deg(thDat(thBothID,3)),rad2deg(thDat(thBothID,6)),binedges);
xlabel('\Delta Theta Angle (Novel - Familiar)')
title("Significant units pre & post");

if saveFlag
    saveas(thetaPie,'subRwdShift_th_pie','png')
    saveas(thMRLFig,'subRwdShift_thMRL_bar','png')
    saveas(thDltMFig,'subRwdShift_thMRL_delta','png')
    saveas(thAngFig,'subRwdShift_thAng_bar','png')
    saveas(thDltAFig,'subRwdShift_thAng_delta','png')
end

%% Theta phase figure
figure; hold on 
plot(rad2deg(0:pi/36:4*pi),cos(pi:pi/36:5*pi),'k','LineWidth',2)
xlim([0 720]); xticks([0 180 360 540 720]); xlabel('Theta Phase')
yticklabels(''); yticks([])
set(gca,'FontSize',12,'FontName','Arial')

%% Spatial Information, peak rate
% lcDat: 1&5 = sig.; 2&6 = SI; 3&7 = pkRate; 4&8 = pkLoc
siFrstID = useCC & lcDat(:,1) <= 0.05;
siLastID = useCC & lcDat(:,5) <= 0.05;
siBothID = siFrstID & siLastID;
nSIBoth = sum(siBothID);

siPie = prepostPie(siFrstID,siLastID,useCC);
title("Sig. Spatial Info. units");

[~,ps.lc_PPSI_both,~,stats.lc_PPSI_both] = ttest(lcDat(siBothID,2), lcDat(siBothID,6));
[~,ps.lc_PPPk_both,~,stats.lc_PPPk_both] = ttest(lcDat(siBothID,3), lcDat(siBothID,7));
[~,ps.lc_PPSI_eith,~,stats.lc_PPSI_eith] = ttest2(lcDat(siFrstID,2),lcDat(siLastID,6));
[~,ps.lc_PPPk_eith,~,stats.lc_PPPk_eith] = ttest2(lcDat(siFrstID,3),lcDat(siLastID,7));

lcBothSIFig = plotBar2(lcDat(siBothID,2),lcDat(siBothID,6));
ylabel('Spatial Information (Bits/spike)')

lcBothPkFig = plotBar2(lcDat(siBothID,3),lcDat(siBothID,7));
ylabel('Peak Field FR (Hz)')

lcEithSIFig = plotBar2(lcDat(siFrstID & ~siLastID,2),lcDat(siLastID & ~siFrstID,6));
ylabel('Spatial Information (Bits/spike)')

lcEithPkFig = plotBar2(lcDat(siFrstID & ~siLastID,3),lcDat(siLastID & ~siFrstID,7));
ylabel('Peak Field FR (Hz)')

if saveFlag
    saveas(siPie,'subRwdShift_SI_pie','png')
    saveas(lcBothSIFig,'subRwdShift_LCSI_both_bar','png')
    saveas(lcBothPkFig,'subRwdShift_LCPk_both_bar','png')
    saveas(lcEithSIFig,'subRwdShift_LCSI_eith_bar','png')
    saveas(lcEithPkFig,'subRwdShift_LCPk_eith_bar','png')
end

%% SI over 10 trial blocks



%% Field peak spatial distribution
binedges = 0:0.05:1.85;
binpos = 0.025:0.05:1.825;

figure; axis square; hold on
plot(100*binpos(lcDat(siBothID,4)),100*binpos(lcDat(siBothID,8)),'ko')
plot([0 100*binpos(end)],[0 100*binpos(end)],'k--')
plot([r1pos r1pos]*100,[0 100*binpos(end)],'r--',[0 100*binpos(end)],[r2pos r2pos]*100,'r--')
xlabel('Absolute Peak Loc. (cm) Familiar')
ylabel('Absolute Peak Loc. (cm) Novel')
set(gca,'FontSize',12,'FontName','Arial')

fieldDstPrePst = histcounts2(binpos(lcDat(siBothID,8)),binpos(lcDat(siBothID,4)),binedges,binedges);
figure; hold on
imagesc(fieldDstPrePst)
colormap("parula")
cbar = colorbar; %clim([prctile(fieldDstPrePst,1,'all'), prctile(fieldDstPrePst,99,'all')]);
plot([0 length(binpos)],[0 length(binpos)],'w--')
plot([r1pos r1pos],[0 binpos(end)],'r--',[0 binpos(end)],[r2pos r2pos],'r--')
xlabel('Absolute Peak Loc. (cm) Familiar'); ylabel('Absolute Peak Loc. (cm) Novel')
ylabel(cbar,'Count','FontSize',12,'Rotation',90)
xticks(1:10:length(binpos))
xticklabels(binpos(1:10:end)*100)
set(gca,'FontSize',12,'FontName','Arial','YDir','normal','XDir','normal')

fieldDstPre = histcounts(binpos(lcDat(siBothID,4)),binedges);
fieldDstPst = histcounts(binpos(lcDat(siBothID,8)),binedges);

prepstFieldDstFig = figure; hold on
plot(100*binpos,fieldDstPre./sum(fieldDstPre),'Color',vColors2(1,:),'LineWidth',2)
plot(100*binpos,fieldDstPst./sum(fieldDstPre),'Color',vColors2(2,:),'LineWidth',2)
plot([r1pos r1pos]*100,[0 0.15],'--','Color',[0.5 0.5 1])
plot([r2pos r2pos]*100,[0 0.15],'--','Color',[0.75 0.75 1])
xlim([0 190]); ylim([0 max([fieldDstPre./sum(fieldDstPre); fieldDstPst./sum(fieldDstPst)+0.02],[],'all')])
xlabel('Track Position (cm)')
ylabel('Probability of field peak')
legend({'Familiar','Novel','RZ1','RZ2'},'Location','northeast')
set(gca,'FontSize',12,'FontName','Arial')

trackLen = max(binedges);
shiftR1 = trackLen/2 - r1pos;
fieldAlignR1 = mod(binpos(lcDat(siBothID,4))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
shiftR2 = trackLen/2 - r2pos;
fieldAlignR2 = mod(binpos(lcDat(siBothID,8))+shiftR2,trackLen)+0.5*dbnsz;
shiftbins = binedges - trackLen/2;

fieldDstPre_Align = histcounts(fieldAlignR1,binedges);
fieldDstPst_Align = histcounts(fieldAlignR2,binedges);

rzDstFig = figure; hold on
plot((binpos-trackLen/2)*100,fieldDstPre_Align./sum(fieldDstPst_Align),'Color',[0.5 0.5 1],'LineWidth',2);
plot((binpos-trackLen/2)*100,fieldDstPst_Align./sum(fieldDstPst_Align),'Color',[0.75 0.75 1],'LineWidth',2);
plot([0 0]*100,[0 0.15],'k--')
ylabel('Probability of field peak')
xlabel('Distance to RZ (cm)')
xlim([-95 95]); ylim([0 max([fieldDstPre_Align./sum(fieldDstPre_Align); fieldDstPst_Align./sum(fieldDstPst_Align)+0.02],[],'all')])
legend({'Familiar','Novel','RZ'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

deltaField_RZAlign = fieldAlignR2 - fieldAlignR1;
circAlignNeg = deltaField_RZAlign < -trackLen/2;
circAlignPos = deltaField_RZAlign > trackLen/2;
deltaField_RZAlign(circAlignNeg) = -(deltaField_RZAlign(circAlignNeg) + trackLen);     % When new field back-shifts
deltaField_RZAlign(circAlignPos) = -(deltaField_RZAlign(circAlignPos) - trackLen);     % When new field forward-shifts

[~,ps.lc_DltLc_both,~,stats.lc_DltLc_both] = ttest(deltaField_RZAlign);

deltaField_Distro = histcounts(deltaField_RZAlign,shiftbins);

fieldShiftRZFig = figure; hold on
bar((shiftbins(1:end-1)+0.5*dbnsz)*100,deltaField_Distro/sum(deltaField_Distro),'FaceColor',[0.25 0.15 1],'HandleVisibility','off')
plot([0 0]*100,[0 max(deltaField_Distro/sum(deltaField_Distro),[],'all')+0.02],'--','Color',[.5 .5 .5],'HandleVisibility','off')
plot(mean(deltaField_RZAlign),max(deltaField_Distro/sum(deltaField_Distro),[],'all')+0.02,'v','Color',[0.25 0.15 1])
xlabel('\Delta RZ-aligned Novel - Familiar (cm)'); ylabel('Probability')
xlim([-trackLen/2-dbnsz trackLen/2+dbnsz]*100)
legend({['mean = ' num2str(mean(deltaField_RZAlign))]})
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(prepstFieldDstFig,'subRwdShift_Field_Distro','png')
    saveas(rzDstFig,'subRwdShift_Field_Align_Distro','png')
    saveas(fieldShiftRZFig,'subRwdShift_Field_Align_Delta','png')
end

%% SPWR Data
swrFrstID = useCC & rpDat(:,2) > 1;
swrLastID = useCC & rpDat(:,4) > 1;
swrBothID = swrFrstID & swrLastID;
nSWRBoth = sum(swrBothID);

tmpSWPie = prepostPie(swrFrstID,swrLastID,useCC);
title("Sig. SPW-R Mod. units");

[~,ps.swr_PPR_both,~,stats.swr_PPR_both] = ttest(rpRat(:,1),rpRat(:,2));
[~,ps.swr_PPP_both,~,stats.swr_PPP_both] = ttest(rpDat(swrBothID,1),rpDat(swrBothID,3));

swPartcpFig = plotBar2(rpDat(swrBothID,1),rpDat(swrBothID,3));
ylabel('P(SPW-R Participation)')

figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
b = bar(mean(rpRat),'FaceColor','flat');
b.CData = vColors2;
errorbar(1:2,mean(rpRat),std(rpRat)/sqrt(size(rpRat,1)),'k.')
plot(rpRat','k-o')
xticks(1:2); xticklabels({'Familiar', 'Novel'})
ylabel('SPW-R Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(tmpSWPie,'subRwdShift_swrMod_pie','png')
    % saveas(siFig,'subRwdShift_SI_bar','png')
    saveas(swPartcpFig,'subRwdShift_swrPartcp_bar','png')
end

%% Functions

function [fhandle] = prepostPie(preSig,postSig,useCCs)

bothSig = preSig & postSig;
notSig = not(preSig | postSig) & useCCs;

nBoth = sum(bothSig);
nPre = sum(preSig);
nPost = sum(postSig);
nNot = sum(notSig);

cMap = [0.25 0.15 1; 0.5 0.5 1; 0.75 0.75 1; 0.25 0.25 0.25];

fhandle = figure;
p = piechart([nBoth, nPre - nBoth, nPost - nBoth, nNot],["Both","Familiar-only","Novel-only","Neither"]);
p.LabelStyle = 'namedata';
colororder(cMap)
end

function [halfStruc] = subEpoch(sess,root,halfStruc)

dbnsz = 0.05;

nLaps = length(sess.valTrials);
epochSubInds = [sess.lapstt(sess.valTrials(1:10:nLaps-mod(nLaps,10))),...
    sess.lapend(sess.valTrials(1:10:nLaps-mod(nLaps,10))+9)];

for i = 1:size(epochSubInds,1)
    [sesstmp, roottmp] = epochStruc(sess,root,epochSubInds(i,:));
    for j = 1:length(root.good)
        cc = root.good(j);
        [halfStruc.subEpochSI(j,i)] = get_SI(roottmp,cc,sesstmp,dbnsz);
    end
end
end

function [fhandle] = plotBar2(dat1,dat2)

vColors2 = [0.5 0.5 1; 0.75 0.75 1];
nUnits = [size(dat1,1) size(dat2,1)];
bardat = [mean(dat1); mean(dat2)];
semdat = [std(dat1)/sqrt(nUnits(1)); std(dat2)/sqrt(nUnits(2))];

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
b = bar([1.15 2.15],bardat,0.3,'FaceColor','flat','BarWidth',0.5);
b.CData = vColors2;
errorbar([1.15 2.15],bardat,semdat,'k.')
violinplot([dat1; dat2],[1*ones(nUnits(1),1); 2*ones(nUnits(2),1)], 'ViolinColor',vColors2,'HalfViolin','left');
xlim([0.5 2.5]); ylim([0 max([dat1; dat2],[],'all')])
xticks(1:2); xticklabels({'Familiar', 'Novel'})
set(gca,'FontSize',12,'FontName','Arial')
end

function [fhandle] = plotDeltaHisto(dat1,dat2,binedges)

deltaDatDistro = histcounts(dat2-dat1,binedges);
bnsz = binedges(2) - binedges(1);

fhandle = figure; hold on
bar(binedges(1:end-1)+bnsz/2,deltaDatDistro/sum(deltaDatDistro),'FaceColor',[0.25 0.15 1],'HandleVisibility','off');
plot(mean(dat2-dat1), max(deltaDatDistro/sum(deltaDatDistro))+0.01,'v','Color',[0.25 0.15 1])
plot([0 0], [0 max(deltaDatDistro/sum(deltaDatDistro))+0.01],'k--')
legend({['mean = ' num2str(mean(dat2-dat1))]})
ylabel('Probability')
set(gca,'FontSize',12,'FontName','Arial')
end