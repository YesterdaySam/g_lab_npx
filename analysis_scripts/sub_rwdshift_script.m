% Subiculum RZ Shift
% spath = 'D:\Data\Kelton\analyses\KW038\KW038_04182025_rec_D5_LLat2';
% spath = 'D:\Data\Kelton\analyses\KW040\KW040_04292025_rec_D4_LLat1';
% spath = 'D:\Data\Kelton\analyses\KW048\KW048_06172025_rec_D2_RLat2';
% spath = 'D:\Data\Kelton\analyses\KW049\KW049_06172025_rec_D2_RLat2';
% spath = 'D:\Data\Kelton\analyses\FG042\FG042_20250607_R3';

% CA1 RZ Shift
% spath = 'D:\Data\Kelton\analyses\FG044\FG044_20250710_R2';

% Sub RZ Shift
% spath = 'D:\Data\Kelton\analyses\KW040\KW040_05012025_rec_D6_LMed1';
% spath = 'D:\Data\Kelton\analyses\KW049\KW049_06202025_rec_D5_LLat2';
% spath = 'D:\Data\Kelton\analyses\KW049\KW049_06222025_rec_D6_LMed1';


cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_RwdShift_Data2.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

rwdShift = sess.valTrials(find(diff(sess.pos(sess.rwdind)) > 0.1,1)+1);   % Find lap of reward shift
nUnits = length(root.good);
saveFlag = 1;

dbnsz = 0.05;
histoBnsz = 5;
wlen = 150;
ripRef = root.ripRef;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
binpos = 0.025:dbnsz:1.825;

% Epoch into pre/post reward shift for reward shift sessions

if ~isempty(rwdShift)
    frstHalfInds         = [sess.ind(1) sess.lapend(rwdShift-1)];
    lastHalfInds         = [sess.lapstt(rwdShift) sess.ind(end)];
    [sessFrst, rootFrst] = epochStruc(sess,root,frstHalfInds);
    [sessLast, rootLast] = epochStruc(sess,root,lastHalfInds);
end

%% Plot example pre/post behavior

lickDIFig = plot_lickDiscrim(sess,[r1pos r2pos]*100,30,10);
plot([rwdShift rwdShift], [-1 1], 'k--','HandleVisibility','off')

bhvrFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,1,6);
ppVelFig = plot_prepost_vel(sessFrst,sessLast);
ppLckFig = plot_prepost_lck(sessFrst,sessLast);

if saveFlag
    saveas(lickDIFig,[root.name '_lickDI'], 'png')
    saveas(bhvrFig,[root.name '_bhvrPrePost'], 'png')
    saveas(ppVelFig,[root.name '_prepostVel'], 'png')
    saveas(ppLckFig,[root.name '_prepostLck'], 'png')
end

%% Plot example pre/post unit
cc = 323;
close all
rzPosFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,1);
rzVelFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,2);
rzThMFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,3);
rzRwdFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,4);
rzSWRFig = plot_prepost(rootFrst,sessFrst,rootLast,sessLast,cc,7);

if saveFlag
    tblind = find(root.info.cluster_id == cc);
    sbase = ['unit_' num2str(cc) '_shank_' num2str(root.info.shankID(tblind)) '_rwdshift_'];
    saveas(rzPosFig,[sbase 'spatial'], 'png')
    saveas(rzVelFig,[sbase 'velocity'],'png')
    saveas(rzThMFig,[sbase 'thetaMod'],'png')
    saveas(rzRwdFig,[sbase 'rewardT'], 'png')
    saveas(rzSWRFig,[sbase 'spwrMod'], 'png')
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
    [frstHalf.trueSI(i),~,frstHalf.truePk(i),frstHalf.trueLc(i),~,~,frstHalf.posfr(i,:),frstHalf.binedges] = get_SI(rootFrst,cc,sessFrst,dbnsz);
    [lastHalf.trueSI(i),~,lastHalf.truePk(i),lastHalf.trueLc(i),~,~,lastHalf.posfr(i,:),lastHalf.binedges] = get_SI(rootLast,cc,sessLast,dbnsz);
    frstHalf.frStandRun(i,:) = get_frStandVRun(rootFrst,cc,sessFrst);
    lastHalf.frStandRun(i,:) = get_frStandVRun(rootLast,cc,sessLast);
    [~,~,frstHalf.trueVelMdl(i)] = plot_frXvel(rootFrst,cc,sessFrst,2,0);
    [~,~,lastHalf.trueVelMdl(i)] = plot_frXvel(rootLast,cc,sessLast,2,0);
    [frstHalf.thetastats(i),frstHalf.thetafr(i,:)] = plot_thetaMod(rootFrst,cc,lfpInd,2*pi/36,0);
    [lastHalf.thetastats(i),lastHalf.thetafr(i,:)] = plot_thetaMod(rootLast,cc,lfpInd,2*pi/36,0);
    frstHalf.swrfr(i,:) = plot_frXripple(rootFrst,cc,sessFrst,ripRef,wlen,histoBnsz,0);
    lastHalf.swrfr(i,:) = plot_frXripple(rootLast,cc,sessLast,ripRef,wlen,histoBnsz,0);
    [~,frstHalf.rwdfr(i,:),frstHalf.trueRI(i)] = plot_frXrwdtime(rootFrst,cc,sessFrst,0.25,5,0);
    [~,lastHalf.rwdfr(i,:),lastHalf.trueRI(i)] = plot_frXrwdtime(rootLast,cc,sessLast,0.25,5,0);
end

frstHalf = subEpochSI(sessFrst,rootFrst,frstHalf,10);
lastHalf = subEpochSI(sessLast,rootLast,lastHalf,10);

frstHalf.ripRate = size(rootFrst.ripStruc(ripRef).ripples,1) / (sum(not(sessFrst.runInds)) / sess.samprate);    %Normalize based on standing periods
lastHalf.ripRate = size(rootLast.ripStruc(ripRef).ripples,1) / (sum(not(sessLast.runInds)) / sess.samprate);
frstHalf.binpos = frstHalf.binedges(1:end-1)+0.5*dbnsz;
lastHalf.binpos = lastHalf.binedges(1:end-1)+0.5*dbnsz;

%% Add lap by lap FR info to epochStrucs

for i = 1:nUnits
    cc = root.good(i);
    [~,frstHalf.frMap(:,:,i)] = get_frXpos(rootFrst,cc,sessFrst,0.05,1.85,1);
    [~,lastHalf.frMap(:,:,i)] = get_frXpos(rootLast,cc,sessLast,0.05,1.85,1);
end

%% Behavioral comparison
% Lick discrimination
[frstHalf.prepLck,frstHalf.lckDI] = get_lickDiscrim(sessFrst,[r1pos r2pos]*100);
[lastHalf.prepLck,lastHalf.lckDI] = get_lickDiscrim(sessLast,[r1pos r2pos]*100);
frstHalf.ulckDI = mean(frstHalf.lckDI,'omitnan');
lastHalf.ulckDI = mean(lastHalf.lckDI,'omitnan');

% Velocity discrimination
frstHalf.preRZV = get_periRZVel(sessFrst,[r1pos r2pos]*100);
lastHalf.preRZV = get_periRZVel(sessLast,[r1pos r2pos]*100);
frstHalf.uPreRZV = mean(frstHalf.preRZV,'omitnan');
lastHalf.uPreRZV = mean(lastHalf.preRZV,'omitnan');

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

        % [~,~,frstHalf.shufRI(i,j)] = plot_frXrwdtime(shiftFrst,cc,shiftSessFrst,0.25,5,0);
        % [~,~,lastHalf.shufRI(i,j)] = plot_frXrwdtime(shiftLast,cc,shiftSessLast,0.25,5,0);
    end
end

toc

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

frStndRunFig = figure; hold on
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

lcBothSIFig = plotBar2(frstHalf.trueSI(siUnits), lastHalf.trueSI(siUnits));
ylabel('Spatial Information (Bits/spike)')

if saveFlag
    sbase = root.name;
    saveas(frStndRunFig,[sbase '_RwdShift_FRStandRun.png'])
    saveas(lcBothSIFig,[sbase '_RwdShift_SI.png'])
end

%% Peak Location distribution
fieldDistroFrst = histcounts(frstHalf.binpos(frstHalf.trueLc(siUnits)),frstHalf.binedges);
fieldDistroLast = histcounts(lastHalf.binpos(lastHalf.trueLc(siUnits)),lastHalf.binedges);

peakDistroFig = plotDistroHisto(fieldDistroFrst,fieldDistroLast,binpos,[r1pos r2pos]);
xlabel('Track Position (cm)')

% Align field distribution around reward
trackLen = max(frstHalf.binedges);
shiftR1 = trackLen/2 - r1pos;
fieldAlignR1 = mod(frstHalf.binpos(frstHalf.trueLc(siUnits))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
shiftR2 = trackLen/2 - r2pos;
fieldAlignR2 = mod(lastHalf.binpos(lastHalf.trueLc(siUnits))+shiftR2,trackLen)+0.5*dbnsz;

fieldDistroFrst_Align = histcounts(fieldAlignR1,frstHalf.binedges);
fieldDistroLast_Align = histcounts(fieldAlignR2,lastHalf.binedges);

rzDistroFig = plotDistroHisto(fieldDistroFrst_Align,fieldDistroLast_Align,binpos-trackLen/2,0);
xlabel('Distance to RZ (cm)')

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
sbase = 'subRwdShift_';

dbnsz = 0.05;
histoBnsz = 5;
binpos = 0.025:dbnsz:1.825;
wlen = 150;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

%%
load('D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\subRwdShift_data2.mat')

nMice = length(siStr);
mID = unique(recID(:,1));
nTotal = size(recID,1);

binedges = 0:5:185;
binpos = 0.025:dbnsz:1.825;
r1posInd = find(binpos > r1pos,1);
r2posInd = find(binpos > r2pos,1);

%% 
frDat = [];     % [frstHalf.standFR, frstHalf.runFR, lastHalf.standFR, lastHalf.runFR]
recID = [];     % [mouseID, recDay, recUnit ID]
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
        elseif subDatT{i,6+j}{1} == 'bdr'
            roiUnits(tmpUnits) = 0;
        end
    end

    pkExclude = frstHalf.truePk < 1 & lastHalf.truePk < 1;
    useUnits = root.info.lyrID(root.goodind) == 1 & root.info.fr(root.goodind) > 0.1 & root.info.uType(root.goodind) & roiUnits(root.goodind) & ~pkExclude;
    nCCs = length(root.good);

    % === Concatenate recording data ===
    useCC = logical([useCC; useUnits]);
    recID = [recID; str2num(subDatT.mouse{i}(end-2:end))*ones(nCCs,1), subDatT.session(i)*ones(nCCs,1), root.good];

    % === Concatenate Velocity-FR data ===
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
    thMap = [thMap; frstHalf.thetafr lastHalf.thetafr];

    % === Concatenate SI and Peak data ===
    lcDat = [lcDat; frstHalf.sig, frstHalf.trueSI, frstHalf.truePk, frstHalf.trueLc,...
        lastHalf.sig, lastHalf.trueSI, lastHalf.truePk, lastHalf.trueLc];
    lcMap = [lcMap; frstHalf.posfr, lastHalf.posfr];

    % === Concatenate SPWR Modulation data ===
    rpDat = [rpDat; frstHalf.ripParticipation', frstHalf.ripModBinCt',...
        lastHalf.ripParticipation', lastHalf.ripModBinCt'];
    rpRat = [rpRat; frstHalf.ripRate, lastHalf.ripRate];
    rpMap = [rpMap; frstHalf.swrfr lastHalf.swrfr];

    % === Concatenate SI over 10-trial blocks ===
    siStr(ct).preBlockSI = frstHalf.subEpochSI;
    siStr(ct).pstBlockSI = lastHalf.subEpochSI;

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

cd(groupSDir)

save([sbase 'data3'],'frDat','recID','useCC','siStr','lcDat','lcMap','thDat','thMap','rpDat','rpRat','rpMap','vlDat','bvDat')

%% Behavior comparisons

uLckDI = [vertcat(bvDat.uPreLckDI), vertcat(bvDat.uPstLckDI)];
uRZVel = [vertcat(bvDat.uPreRZVel), vertcat(bvDat.uPstRZVel)];

[~,ps.bhv_PPLckDI,~,stats.bhv_PPLckDI] = ttest(uLckDI(:,1),uLckDI(:,2));
[~,ps.bhv_PPRZ1,~,stats.bhv_PPRZ1] = ttest(uRZVel(:,1),uRZVel(:,3));
[~,ps.bhv_PPRZ2,~,stats.bhv_PPRZ2] = ttest(uRZVel(:,2),uRZVel(:,4));

uLckDIFig = plotBarByMouse(uLckDI);
ylabel('Lick DI ((RZ - AZ) / (RZ + AZ))'); ylim([-1 1])

uRZ1VlFig = plotBarByMouse(uRZVel(:,[1,3]));
ylabel('Velocity (cm/s) 30cm Pre RZ1')

uRZ2VlFig = plotBarByMouse(uRZVel(:,[2,4]));
ylabel('Velocity (cm/s) 30cm Pre RZ2')

if saveFlag
    saveas(uLckDIFig,[sbase 'bhv_LckDI_bar'],'png')
    saveas(uRZ1VlFig,[sbase 'bhv_RZ1Vl_bar'],'png')
    saveas(uRZ2VlFig,[sbase 'bhv_RZ2Vl_bar'],'png')
end

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
vlDltBFig = plotDeltaHisto(vlDat(vlBothID,2),vlDat(vlBothID,5),binedges);
xlabel('\Delta Slope (Novel - Familiar)')
title("Significant units pre & post");

vlDltRFig = plotDeltaHisto(vlDat(vlBothID,3),vlDat(vlBothID,6),binedges);
ylabel('Probability'); xlabel('\Delta R (Novel - Familiar)')
title("Significant units pre & post");

vlBothCounts = [groupcounts(recID(vlFrstID,1)) groupcounts(recID(vlLastID,1))];
[~,ps.vl_ModCt_both,~,stats.vl_ModCt_both] = ttest(vlBothCounts(:,1),vlBothCounts(:,2));
vlBothCtFig = plotBarByMouse(vlBothCounts);
ylabel('# Sig. Velocity')

if saveFlag
    saveas(velPieFig,[sbase 'vel_ModCt_pie'],'png')
    saveas(vlDltBFig,[sbase 'vel_deltaB'],'png')
    saveas(vlDltRFig,[sbase 'vel_deltaR'],'png')
    saveas(vlBothCtFig,[sbase 'vel_ModCt_bar'],'png')
end

%% Firing Rate stand vs run
% frDat: 1&3 = standing; 2&4 = running
nUseCC = sum(useCC);

[~,ps.fr_PPStnd_all,~,stats.fr_PPStnd_all] = ttest(frDat(useCC,1),frDat(useCC,3));
[~,ps.fr_PPRunn_all,~,stats.fr_PPRunn_all] = ttest(frDat(useCC,2),frDat(useCC,4));

bardat = [mean(frDat(useCC,1)), mean(frDat(useCC,2)); mean(frDat(useCC,3)), mean(frDat(useCC,4))];
xcoords = ones(sum(useCC),1);

frStndRunFig = figure; hold on
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
    saveas(frStndRunFig,[sbase 'fr_bar'],'png')
end

%% Theta Modulation
% thDat: 1&4 = sig.; 2&5 = MRL; 3&6 = Angle
thFrstID = useCC & thDat(:,1) <= 0.05;
thLastID = useCC & thDat(:,4) <= 0.05;
thBothID = thFrstID & thLastID;

thPieFig = prepostPie(thFrstID,thLastID,useCC);
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

thBothCounts = [groupcounts(recID(thFrstID,1)) groupcounts(recID(thLastID,1))];
[~,ps.th_ModCt_both,~,stats.th_ModCt_both] = ttest(thBothCounts(:,1),thBothCounts(:,2));
thBothCtFig = plotBarByMouse(thBothCounts);
ylabel('# Sig. Theta-Mod')

if saveFlag
    saveas(thPieFig,[sbase 'th_Mod_pie'],'png')
    saveas(thMRLFig,[sbase 'th_MRL_bar'],'png')
    saveas(thDltMFig,[sbase 'th_MRL_delta'],'png')
    saveas(thAngFig,[sbase 'th_Ang_bar'],'png')
    saveas(thDltAFig,[sbase 'th_Ang_delta'],'png')
    saveas(thBothCtFig,[sbase 'th_ModCt_bar'],'png')
end

%% Waterfall by theta phase
binedges = rad2deg(0:pi/18:2*pi);
thPeak = find(binedges == 180,1);

[thBothPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(thMap(thBothID,1:length(binedges)-1),binedges);
plot([thPeak thPeak],[0 sum(thBothID)],'w--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Theta Phase')
thBothPreSortPreHisto = plot_unitPkHisto(tmpMap,binedges);
xlabel('Theta Phase')

[thBothPstSortPstFig,tmpMap] = plot_unitWaterfall(thMap(thBothID,length(binedges):end),binedges);
plot([thPeak thPeak],[0 sum(thBothID)],'w--','LineWidth',2)
title('Novel RZ, sort Novel'); xlabel('Theta Phase')
thBothPstSortPstHisto = plot_unitPkHisto(tmpMap,binedges);
xlabel('Theta Phase')

[thBothPstSortPreFig] = plot_unitWaterfall(thMap(thBothID,length(binedges):end),binedges,sortPre);
plot([thPeak thPeak],[0 sum(thBothID)],'w--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Theta Phase')

if saveFlag
    saveas(thBothPreSortPreFig,[sbase 'th_both_pre_sortPre'],'png')
    saveas(thBothPstSortPstFig,[sbase 'th_both_pst_sortPst'],'png')
    saveas(thBothPstSortPreFig,[sbase 'th_both_pst_sortPre'],'png')
    saveas(thBothPreSortPreHisto,[sbase 'th_both_pre_distro'],'png')
    saveas(thBothPstSortPstHisto,[sbase 'th_both_pst_distro'],'png')
end

%% Theta phase figure
figure; hold on
cycleMax = 2*pi;
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.12])
% plot(rad2deg(0:pi/36:cycleMax),cos(pi:pi/36:5*pi),'k','LineWidth',2)    % Trough = 0/360
plot(rad2deg(0:pi/36:cycleMax),cos(0:pi/36:cycleMax),'k','LineWidth',2)     % Trough = 180
xlim([0 rad2deg(cycleMax)]); xticks([0 180 360 540 720]); xlabel('Theta Phase')
yticklabels(''); yticks([])
set(gca,'FontSize',12,'FontName','Arial')

%% Spatial Information, peak rate
% lcDat: 1&5 = sig.; 2&6 = SI; 3&7 = pkRate; 4&8 = pkLoc
siFrstID = useCC & lcDat(:,1) <= 0.05;
siLastID = useCC & lcDat(:,5) <= 0.05;
siBothID = siFrstID & siLastID;

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

siBothCounts = [groupcounts(recID(siFrstID,1)) groupcounts(recID(siLastID,1))];
siBothRatio = siBothCounts ./ [groupcounts(recID(~siFrstID,1)) groupcounts(recID(~siLastID,1))];
[~,ps.lc_SICt_both,~,stats.lc_SICt_both] = ttest(siBothCounts(:,1),siBothCounts(:,2));
lcBothCtFig = plotBarByMouse(siBothRatio);
ylabel('P(Units with significant Spatial Info.)')

if saveFlag
    saveas(siPie,[sbase 'si_Mod_pie'],'png')
    saveas(lcBothSIFig,[sbase 'si_both_bar'],'png')
    saveas(lcBothPkFig,[sbase 'pk_both_bar'],'png')
    saveas(lcEithSIFig,[sbase 'si_eith_bar'],'png')
    saveas(lcEithPkFig,[sbase 'pk_eith_bar'],'png')
    saveas(lcBothCtFig,[sbase 'si_Mod_both_bar'],'png')
end

%% Waterfall by position

[spBothPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(lcMap(siBothID,1:length(binpos)),binedges);
plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
spBothPreSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])

[spBothPstSortPstFig,tmpMap] = plot_unitWaterfall(lcMap(siBothID,length(binpos)+1:end),binedges);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Novel'); xlabel('Track Position (cm)')
spBothPstSortPstHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])
%%
[spBothPstSortPreFig,tmpMap] = plot_unitWaterfall(lcMap(siBothID,length(binpos)+1:end),binedges,sortPre);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')
spBothPstSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])
%%

[spEithPreSortPreFig,~,sortPre] = plot_unitWaterfall(lcMap(siFrstID & ~siLastID,1:length(binpos)),binedges);
plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')

[spEithPstSortPreFig] = plot_unitWaterfall(lcMap(siFrstID & ~siLastID,length(binpos)+1:end),binedges,sortPre);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')

[spEithPstSortPstFig,~,Sortpst] = plot_unitWaterfall(lcMap(siLastID & ~siFrstID,length(binpos)+1:end),binedges);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Novel'); xlabel('Track Position (cm)')

[spEithPreSortPstFig] = plot_unitWaterfall(lcMap(siLastID & ~siFrstID,1:length(binpos)),binedges,Sortpst);
plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Novel'); xlabel('Track Position (cm)')

if saveFlag
    saveas(spBothPreSortPreFig,[sbase 'sp_both_pre_sortPre'],'png')
    saveas(spBothPstSortPstFig,[sbase 'sp_both_pst_sortPst'],'png')
    saveas(spBothPstSortPreFig,[sbase 'sp_both_pst_sortPre'],'png')
    saveas(spEithPreSortPreFig,[sbase 'sp_eith_pre_sortPre'],'png')
    saveas(spEithPstSortPreFig,[sbase 'sp_eith_pst_sortPre'],'png')
    saveas(spEithPstSortPstFig,[sbase 'sp_eith_pst_sortPst'],'png')
    saveas(spEithPreSortPstFig,[sbase 'sp_eith_pre_sortPst'],'png')
end

%% SI over 10 trial blocks

cmapcool = cool(nMice);
siEpochFig = figure; hold on
for i = 1:nMice
    subID = siBothID(recID(:,1) == mID(i));
    siStr(i).meanPreSI = mean(siStr(i).preBlockSI(subID,:));
    nEpochsPre(i) = size(siStr(i).preBlockSI,2);
    plot(1:nEpochsPre(i),siStr(i).meanPreSI,'-o','Color',cmapcool(i,:))
    nEpochsPst(i) = size(siStr(i).pstBlockSI,2);
end
maxEpoch = max(nEpochsPre);
uPreMat = NaN(size(recID,1),maxEpoch);
uPstMat = NaN(size(recID,1),max(nEpochsPst));

ct = 1;
for i = 1:nMice
    nUnits = sum(recID(:,1) == mID(i));
    uPreMat(ct:nUnits+ct-1,1:nEpochsPre(i)) = siStr(i).preBlockSI;
    uPstMat(ct:nUnits+ct-1,1:nEpochsPst(i)) = siStr(i).pstBlockSI;
    ct = ct+nUnits;
end
for i = 1:nMice
    subID = siBothID(recID(:,1) == mID(i));
    siStr(i).meanPstSI = mean(siStr(i).pstBlockSI(subID,:));
    nPst(i) = size(siStr(i).pstBlockSI,2);
    plot(maxEpoch+1:nPst(i)+maxEpoch,siStr(i).meanPstSI,'-o','Color',cmapcool(i,:))
end

errorbar(1:maxEpoch,mean(uPreMat(siBothID,:),'omitnan'),std(uPreMat(siBothID,:),'omitnan')./sqrt(nTotal),'k-','LineWidth',2)
errorbar(maxEpoch+1:max(nEpochsPst)+maxEpoch,mean(uPstMat(siBothID,:),'omitnan'),std(uPstMat(siBothID,:),'omitnan')./sqrt(nTotal),'k-','LineWidth',2)
plot([maxEpoch+0.5 maxEpoch+0.5],[0 max(ylim)],'k--')
xlabel('10-Trial block #')
ylabel('Spatial Information (bits/spike)')
% ylim([0 1])
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(siEpochFig,[sbase 'si_both_subEpoch'],'png')
end

%% Field peak spatial distribution

% Familiar vs Novel RZ correlation
xrand = 2*rand(size(lcDat(siBothID,4)))-1;
yrand = 2*rand(size(lcDat(siBothID,4)))-1;
pkPosFig = figure; axis square; hold on
plot(100*binpos(lcDat(siBothID,4))+xrand',100*binpos(lcDat(siBothID,8))+yrand','k.','MarkerSize',10)
plot([0 100*binpos(end)],[0 100*binpos(end)],'k--')
plot([r1pos r1pos]*100,[0 100*binpos(end)],'r--',[0 100*binpos(end)],[r2pos r2pos]*100,'r--')
xlabel('Absolute Peak Loc. (cm) Familiar'); xlim([0 100*binpos(end)])
ylabel('Absolute Peak Loc. (cm) Novel');    ylim([0 100*binpos(end)])
set(gca,'FontSize',12,'FontName','Arial')

% Familiar vs Novel RZ correlation heatmap
fieldDstPrePst = histcounts2(binpos(lcDat(siBothID,8)),binpos(lcDat(siBothID,4)),binedges,binedges);
pkPosMapFig = figure; hold on; axis square
imagesc(fieldDstPrePst)
colormap("parula")
cbar = colorbar; %clim([prctile(fieldDstPrePst,1,'all'), prctile(fieldDstPrePst,99,'all')]);
plot([0 length(binpos)],[0 length(binpos)],'w--')
plot([r1posInd r1posInd],[0 length(binpos)],'r--',[0 length(binpos)],[r2posInd r2posInd],'r--')
xlabel('Absolute Peak Loc. (cm) Familiar'); ylabel('Absolute Peak Loc. (cm) Novel')
ylabel(cbar,'Count','FontSize',12,'Rotation',90); 
xlim([0 length(binedges)]); ylim([0 length(binedges)])
xticks(1:10:length(binpos)); yticks(1:10:length(binpos))
xticklabels(binedges(1:10:end)*100); yticklabels(binedges(1:10:end)*100)
set(gca,'FontSize',12,'FontName','Arial','YDir','normal','XDir','normal')

%% Finding Track Relative or Reward Relative cells
spPk1 = binpos(lcDat(:,4));
spPk2 = binpos(lcDat(:,8));

dth = 0.3;  % Distance from peak threshold (m)
drz = r2pos - r1pos;
trCells = (spPk2 - spPk1 < dth & spPk2 - spPk1 > -dth)' & siBothID; % threshold 30cm
rrCells = ((spPk2 - spPk1 > drz-dth & spPk2 - spPk1 < drz+dth)' & siBothID) | ((spPk2 - spPk1 < -drz+dth & spPk2 - spPk1 > -drz-dth)' & siBothID);
irCells = siBothID & ~trCells & ~rrCells;

figure; hold on; axis square
patch(100*[0 dth 1.85 1.85 1.85-dth 0 0],100*[0 0 1.85-dth 1.85 1.85 dth 0],...
    'b','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
patch(100*[0 1.25 0.65 0 0],100*[0.6 1.85 1.85 1.2 0.6],...
    'r','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
patch(100*[1.2 1.85 1.85 0.6 1.2],100*[0 0.65 1.25 0 0],...
    'r','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
plot(100*spPk1(siBothID)+xrand',100*spPk2(siBothID)+yrand','k.','MarkerSize',10)
% plot(100*spPk1(irCells),100*spPk2(irCells),'k.')
% plot(100*spPk1(trCells),100*spPk2(trCells),'b.')
% plot(100*spPk1(rrCells),100*spPk2(rrCells),'r.')
plot([0 100*binpos(end)],[0 100*binpos(end)],'k--')
plot([r1pos r1pos]*100,[0 100*binpos(end)],'r--',[0 100*binpos(end)],[r2pos r2pos]*100,'r--')
xlabel('Absolute Peak Loc. (cm) Familiar'); xlim([0 100*binpos(end)])
ylabel('Absolute Peak Loc. (cm) Novel');    ylim([0 100*binpos(end)])
set(gca,'FontSize',12,'FontName','Arial')

%% Distributions over linearized space
fieldDstPre = histcounts(spPk1,binedges);
fieldDstPst = histcounts(spPk2,binedges);

prepstFieldDstFig = plotDistroHisto(fieldDstPre,fieldDstPst,binpos,[r1pos r2pos]);
xlabel('Track Position (cm)')

trackLen = max(binedges);
shiftR1 = trackLen/2 - r1pos;
fieldAlignR1 = mod(binpos(lcDat(siBothID,4))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
shiftR2 = trackLen/2 - r2pos;
fieldAlignR2 = mod(binpos(lcDat(siBothID,8))+shiftR2,trackLen)+0.5*dbnsz;
shiftbins = binedges - trackLen/2;

fieldDstPre_Align = histcounts(fieldAlignR1,binedges);
fieldDstPst_Align = histcounts(fieldAlignR2,binedges);

rzDstFig = plotDistroHisto(fieldDstPre_Align,fieldDstPst_Align,binpos-trackLen/2,0);
xlabel('Distance to RZ (cm)')

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

% Shuffle novel RZ peak locations and calculate global confidence bands
clear deltaField_DistroJit
for i = 1:250
    fieldAlignR1 = mod(binpos(lcDat(siBothID,4))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
    % fieldAbsR1 = binpos(lcDat(siBothID,4));

    rShift = 1 + randi(length(binedges) - 2,sum(siBothID),1);
    jitPost = mod(lcDat(siBothID,8) + rShift, length(binedges)-1)+1;
    fieldAlignJit = mod(binpos(jitPost)+shiftR2,trackLen)+0.5*dbnsz;
    fieldDstJit_Align = histcounts(fieldAlignJit,binedges);
    deltaField_RZJit = fieldAlignJit - fieldAlignR1;
    circAlignNeg = deltaField_RZJit < -trackLen/2;
    circAlignPos = deltaField_RZJit > trackLen/2;
    deltaField_RZJit(circAlignNeg) = -(deltaField_RZJit(circAlignNeg) + trackLen);     % When new field back-shifts
    deltaField_RZJit(circAlignPos) = -(deltaField_RZJit(circAlignPos) - trackLen);     % When new field forward-shifts
    deltaField_DistroJit(:,i) = histcounts(deltaField_RZJit,shiftbins);
    % absField_DistroJit(:,i) = histcounts(binpos(jitPost)-fieldAbsR1,shiftbins);
end
[~,fieldShiftRZJitFig] = get_confband(deltaField_DistroJit',deltaField_Distro,1,shiftbins(1:end-1)*100,dbnsz*100);
xlabel('\Delta RZ-aligned Novel - Familiar (cm)');

if saveFlag
    saveas(pkPosFig,[sbase 'lc_PeakComp'],'png')
    saveas(pkPosMapFig,[sbase 'lc_PeakMap'],'png')
    saveas(prepstFieldDstFig,[sbase 'lc_Distro'],'png')
    saveas(rzDstFig,[sbase 'lc_Align_Distro'],'png')
    saveas(fieldShiftRZFig,[sbase 'lc_Align_Delta'],'png')
    saveas(fieldShiftRZJitFig,[sbase 'lc_Align_Delta_shuf'],'png')
end

%% Waterfall by TR or RR

cMap = [0.25 0.15 1; 0.75 0.75 1; 0.25 0.25 0.25];

spTrRrPie = figure;
p = piechart([sum(trCells), sum(rrCells), sum(siBothID)-sum(trCells)-sum(rrCells)],["Track-Relative","Reward-Relative","Intermediate"]);
p.LabelStyle = 'namepercent';
colororder(cMap)

[spTRPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(lcMap(trCells,1:length(binpos)),binedges);
plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
spTRPreSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])

[spTRPstSortPreFig,tmpMap] = plot_unitWaterfall(lcMap(trCells,length(binpos)+1:end),binedges,sortPre);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')
spTRPstSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])

[spRRPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(lcMap(rrCells,1:length(binpos)),binedges);
plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
spRRPreSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])

[spRRPstSortPreFig,tmpMap] = plot_unitWaterfall(lcMap(rrCells,length(binpos)+1:end),binedges,sortPre);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')
spRRPstSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100);
xlabel('Track Position (cm)'); ylim([0 0.15])

% [spIRPreSortPreFig,~,sortPre] = plot_unitWaterfall(lcMap(irCells,1:length(binpos)),binedges);
% plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
% 
% spIRPstSortPreFig = plot_unitWaterfall(lcMap(irCells,length(binpos)+1:end),binedges,sortPre);
% plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')

if saveFlag
    saveas(spTrRrPie,[sbase,'sp_trrr_pie'],'png')
    saveas(spTRPreSortPreFig,[sbase 'sp_tr_pre_sortPre'],'png')
    saveas(spTRPstSortPreFig,[sbase 'sp_tr_pst_sortPre'],'png')
    saveas(spRRPreSortPreFig,[sbase 'sp_rr_pre_sortPre'],'png')
    saveas(spRRPstSortPreFig,[sbase 'sp_rr_pst_sortPre'],'png')
    saveas(spTRPreSortPreHisto,[sbase 'sp_tr_pre_sortPre_histo'],'png')
    saveas(spTRPstSortPreHisto,[sbase 'sp_tr_pst_sortPre_histo'],'png')
    saveas(spRRPreSortPreHisto,[sbase 'sp_rr_pre_sortPre_histo'],'png')
    saveas(spRRPstSortPreHisto,[sbase 'sp_rr_pst_sortPre_histo'],'png')
    % saveas(spIRPreSortPreFig,[sbase 'sp_ir_pre_sortPre'],'png')
    % saveas(spIRPstSortPreFig,[sbase 'sp_ir_pst_sortPre'],'png')
end

%% Population Vector analysis
% Sparseness - P(alpha) From Battaglia et al., 2004
% Normalize all PF firing rates, sum by bin, and sqaure
% Divide by sum of squared PF firing rates by bin, and normalize by N cells

activeUnits = siBothID;
nBins = length(binpos);

posNormPre = lcMap(activeUnits,1:nBins) ./ max(lcMap(activeUnits,1:nBins),[],2);
posNormPst = lcMap(activeUnits,nBins+1:end) ./ max(lcMap(activeUnits,nBins+1:end),[],2);

% numerPre = sum(posNormPre).^2;
% denomPre = sum(posNormPre.^2);
% 
% pvSparsePre = numerPre ./ denomPre ./ size(posNormPre,1);

pvPrePst = corr(posNormPst,posNormPre);

pvFig = figure; hold on; axis square
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
imagesc(pvPrePst,[prctile(pvPrePst,0,'all'), prctile(pvPrePst,99,'all')])
colormap("turbo")
cbarLims = [0 1];
cbar = colorbar('Ticks',cbarLims,'TickLabels',cbarLims); clim(cbarLims); % ,'position', [0.9 0.12 0.05 0.3]
xlim([-0 nBins+0.5]); xticks([1 nBins]); xticklabels({'0','185'}); xlabel('Position (cm) Familiar RZ')
ylim([-0 nBins+0.5]); yticks([1 nBins]); yticklabels({'0','185'}); ylabel('Position (cm) Novel RZ')
% plot([r1posInd r1posInd],[0 nBins],'w--',[0 nBins],[r2posInd r2posInd],'w--')
plot([0 nBins],[0 nBins],'w--');
plot([0 nBins/2],[nBins/2 nBins],'w--',[nBins/2 nBins],[0 nBins/2],'w--');
ylabel(cbar,'Pearson Correlation of PV','FontSize',12,'Rotation',90)
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

% Compare identity to off-diagonal
idMat = logical(eye(size(pvPrePst)));
dgMat = logical(spdiags([1 1],[-round(nBins/2) round(nBins/2)],nBins,nBins));

uPVCorrID = [];
uPVCorrDG = [];

for i = 1:nMice
    tmpUnits = activeUnits & recID(:,1) == mID(i);

    posNormPre = lcMap(tmpUnits,1:nBins) ./ max(lcMap(tmpUnits,1:nBins),[],2);
    posNormPst = lcMap(tmpUnits,nBins+1:end) ./ max(lcMap(tmpUnits,nBins+1:end),[],2);
    try
        pvPrePst = corr(posNormPst,posNormPre);
    catch
        pvPrePst = NaN(nBins,nBins);
    end
    uPVCorrID(i) = mean(pvPrePst(idMat),'all');
    uPVCorrDG(i) = mean(pvPrePst(dgMat),'all');
end

[~,ps.lc_pv_tr_idVdg,~,stats.lc_pv_tr_idVdg] = ttest(uPVCorrID,uPVCorrDG);
ps.lc_pv_tr_idVdg

pvCompFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.4])
plot([uPVCorrID' uPVCorrDG']','-o','Color',[.5 .5 .5])
errorbar([1 2],mean([uPVCorrID' uPVCorrDG']),std([uPVCorrID' uPVCorrDG'])./sqrt(nMice),'k.','LineWidth',2,'CapSize',20)
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Diagonal','Off-Diag'})
ylim([-0.25 1]); ylabel('Mean of PV Correlation')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(pvFig,[sbase 'lc_pv_corr_trCells'],'png')
    saveas(pvCompFig,[sbase 'lc_pv_comp_trCells'],'png')
end

%% PV analysis through time

lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = lyrUnits & hiFRUnits & root.info.uType(root.goodind);
bothSIUnits = useUnits & (lastHalf.sig <= 0.05 & frstHalf.sig <= 0.05);

posNormPre = frstHalf.posfr(bothSIUnits,:) ./ max(frstHalf.posfr(bothSIUnits,:),[],2);
posNormPst = lastHalf.posfr(bothSIUnits,:) ./ max(lastHalf.posfr(bothSIUnits,:),[],2);

for i = 1:size(frstHalf.frMap,1)
    posNormPreLap = squeeze(frstHalf.frMap(i,:,bothSIUnits)) ./ max(squeeze(frstHalf.frMap(i,:,bothSIUnits)),[],1);
    pvPrePre = corr(posNormPreLap',posNormPre,'rows','complete');
    pvPrePst = corr(posNormPreLap',posNormPst,'rows','complete');
    uPVPrePre(i) = mean(pvPrePre(idMat),'all');
    uPVPrePst(i) = mean(pvPrePst(idMat),'all');
end
for i = 1:size(lastHalf.frMap,1)
    posNormPstLap = squeeze(lastHalf.frMap(i,:,bothSIUnits)) ./ max(squeeze(lastHalf.frMap(i,:,bothSIUnits)),[],1);
    pvPstPre = corr(posNormPstLap',posNormPre,'rows','complete');
    pvPstPst = corr(posNormPstLap',posNormPst,'rows','complete');
    uPVPstPre(i) = mean(pvPstPre(idMat),'all');
    uPVPstPst(i) = mean(pvPstPst(idMat),'all');
end

figure; hold on;
plot(1:size(frstHalf.frMap,1),smoothdata(uPVPrePre),'b',1:size(lastHalf.frMap,1),smoothdata(uPVPstPre),'r')
plot(1:size(frstHalf.frMap,1),smoothdata(uPVPrePst),'Color',[.6 .6 1])
plot(1:size(lastHalf.frMap,1),smoothdata(uPVPstPst),'Color',[1 .6 .6]);
ylim([0 1]); xlabel('Lap'); ylabel('Mean PV Corr. to Familiar RZ PV')
legend('Fam->Fam','Novel->Fam','Fam->Nov','Nov->Nov')
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

% pvFig = figure; hold on; axis square
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
% imagesc(pvPreLap,[prctile(pvPreLap,0,'all'), prctile(pvPreLap,99,'all')])
% colormap("turbo")
% cbarLims = [0 1];
% cbar = colorbar('Ticks',cbarLims,'TickLabels',cbarLims); clim(cbarLims); % ,'position', [0.9 0.12 0.05 0.3]
% xlim([-0 nBins+0.5]); xticks([1 nBins]); xticklabels({'0','185'}); xlabel('Position (cm) Familiar RZ')
% ylim([-0 nBins+0.5]); yticks([1 nBins]); yticklabels({'0','185'}); ylabel('Position (cm) Novel RZ')
% % plot([r1posInd r1posInd],[0 nBins],'w--',[0 nBins],[r2posInd r2posInd],'w--')
% plot([0 nBins],[0 nBins],'w--');
% plot([0 nBins/2],[nBins/2 nBins],'w--',[nBins/2 nBins],[0 nBins/2],'w--');
% ylabel(cbar,'Pearson Correlation of PV','FontSize',12,'Rotation',90)
% set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

%% SPWR Data
swrFrstID = useCC & rpDat(:,2) > 1;
swrLastID = useCC & rpDat(:,4) > 1;
swrBothID = swrFrstID & swrLastID;

tmpSWPie = prepostPie(swrFrstID,swrLastID,useCC);
title("Sig. SPW-R Mod. units");

[~,ps.swr_PPR_both,~,stats.swr_PPR_both] = ttest(rpRat(:,1),rpRat(:,2));
[~,ps.swr_PPP_both,~,stats.swr_PPP_both] = ttest(rpDat(swrBothID,1),rpDat(swrBothID,3));

rpRatFig = plotBarByMouse(rpRat);
ylabel('SPW-R Rate (Hz)')

swPartcpFig = plotBar2(rpDat(swrBothID,1),rpDat(swrBothID,3));
ylabel('P(SPW-R Participation)')

swrBothCounts = [groupcounts(recID(swrFrstID,1)) groupcounts(recID(swrLastID,1))];
[~,ps.swr_ModCt_both,~,stats.swr_ModCt_both] = ttest(swrBothCounts(:,1),swrBothCounts(:,2));
swrModCtFig = plotBarByMouse(swrBothCounts);
ylabel('# SPW-R Modulated')

if saveFlag
    saveas(tmpSWPie,[sbase 'swr_Mod_pie'],'png')
    saveas(rpRatFig,[sbase 'swr_Rate_bar'],'png')
    saveas(swPartcpFig,[sbase 'swr_Partcp_bar'],'png')
    saveas(swrModCtFig,[sbase 'swr_ModCt_bar'],'png')
end

%% Waterfall by sharp wave ripple peak
binedges = -wlen:histoBnsz:wlen;
swrPeak = find(binedges == 0,1);

[swrBothPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(rpMap(swrBothID,1:length(binedges)-1),binedges);
plot([swrPeak swrPeak],[0 sum(swrBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Time to SWR peak (ms)')
swrBothPreSortPreHisto = plot_unitPkHisto(tmpMap,binedges);
xlabel('Time to SWR peak (ms)')

[swrBothPstSortPstFig,tmpMap] = plot_unitWaterfall(rpMap(swrBothID,length(binedges):end),binedges);
plot([swrPeak swrPeak],[0 sum(swrBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Novel'); xlabel('Time to SWR peak (ms)')
swrBothPstSortPstHisto = plot_unitPkHisto(tmpMap,binedges);
xlabel('Time to SWR peak (ms)')

[swrBothPstSortPreFig] = plot_unitWaterfall(rpMap(swrBothID,length(binedges):end),binedges,sortPre);
plot([swrPeak swrPeak],[0 sum(swrBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Time to SWR peak (ms)')

if saveFlag
    saveas(swrBothPreSortPreFig,[sbase 'swr_both_pre_sortPre'],'png')
    saveas(swrBothPstSortPstFig,[sbase 'swr_both_pst_sortPst'],'png')
    saveas(swrBothPstSortPreFig,[sbase 'swr_both_pst_sortPre'],'png')
    saveas(swrBothPreSortPreHisto,[sbase 'swr_both_pre_distro'],'png')
    saveas(swrBothPstSortPstHisto,[sbase 'swr_both_pst_distro'],'png')
end

%% Save stats

if saveFlag
    save([sbase, 'stats'],'ps','stats')
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

function [fhandle] = plotBarByMouse(dat)
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.12 0.39])
b = bar(mean(dat),'FaceColor','flat','HandleVisibility','off');
b.CData = vColors2;
errorbar(1:2,mean(dat),std(dat)/sqrt(size(dat,1)),'k.','HandleVisibility','off')
plot(dat','k-o')
xlim([0.5 2.5])
xticks(1:2); xticklabels({'Familiar', 'Novel'})
legend('Mouse')
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

function [fhandle] = plotDistroHisto(distro1,distro2,binpos,rzPos)
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

fhandle = figure; hold on
plot(binpos*100,distro1./sum(distro1),'Color',vColors2(1,:),'LineWidth',2);
plot(binpos*100,distro2./sum(distro2),'Color',vColors2(2,:),'LineWidth',2);
if length(rzPos) == 2
    plot([rzPos(1) rzPos(1)]*100,[0 0.15],'--','Color',vColors2(1,:))
    plot([rzPos(2) rzPos(2)]*100,[0 0.15],'--','Color',vColors2(2,:))
else
    plot([rzPos rzPos]*100,[0 0.15],'k--')
end
ylim([0 max([distro1./sum(distro1); distro2./sum(distro2)+0.02],[],'all')])
xlim([binpos(1)*100-1 binpos(end)*100+1])
ylabel('Probability of field peak')
legend({'Familiar','Novel'},'Location','northeast')
set(gca,'FontSize',12,'FontName','Arial')
end

function [fhandle] = plot_unitPkHisto(frMapRaw,binedges)
nBins = size(frMapRaw,2);
nUnits = size(frMapRaw,1);

unitMax = max(frMapRaw,[],2);
frMapPk = zeros(size(frMapRaw));

for i = 1:nUnits
    tmpbns = find(frMapRaw(i,:) == unitMax(i)); % In case of multiple peak normalized bins
    frMapPk(i,tmpbns(1)) = 1;
end

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
plot(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(frMapPk)./sum(frMapPk,'all'),'Color',[0.25 0.15 1])
% bar(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(frMapPk)./sum(frMapPk,'all'),'FaceColor',[0.25 0.15 1])
xlim([0 max(binedges)])
ylabel('Probability')
set(gca,'FontSize',12,'FontName','Arial')
end