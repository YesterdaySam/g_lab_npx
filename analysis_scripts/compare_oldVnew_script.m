%% Prepare new data
clear
close all

% spath = 'D:\Data\Kelton\analyses\ZM002\ZM002_09302025_rec_D4_LLat1';
% spath = 'D:\Data\Kelton\analyses\KW073\KW073_10162026_rec_D3_LLat1';
spath = 'D:\Data\Kelton\analyses\ZM006\ZM006_10162026_rec_D3_LLat1';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
datfile = dir("*_datStruc.mat");
try load(datfile.name); catch; disp('No existing shuffled data file'); end

nUnits = length(root.good);
saveFlag = 1;

dbnsz = 0.05;
histoBnsz = 5;
vspbnsz = 0.01;
wlen = 150;
ripRef = root.ripRef;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
binpos = 0.025:dbnsz:1.825;
sbase = root.name;

%% Arrange behavior

clear datStruc
datStruc.name   = root.name;

lickDIFig = plot_lickDiscrim(sess,[r1pos r2pos]*100,15,10);

[~,~,vXCorFig] = plot_lap_velCCorr(sess);
datStruc.spVelMap = plot_lap_velCCorr(sess,vspbnsz,0.002,0);

datStruc.vCorXLap = corr(datStruc.spVelMap',mean(datStruc.spVelMap)','rows','complete');
tmp = [datStruc.vCorXLap; NaN(5-mod(size(datStruc.vCorXLap,1),5),1)]; 
datStruc.vCorAvg = nanmean(reshape(tmp,5,size(tmp,1)/5));

% Lick discrimination
[datStruc.prepLck,datStruc.lckDI] = get_lickDiscrim(sess,[r1pos r2pos]*100);
datStruc.ulckDI = mean(datStruc.lckDI,'omitnan');

% Velocity discrimination
datStruc.preRZV = get_periRZVel(sess,[r1pos r2pos]*100);
datStruc.uPreRZV = mean(datStruc.preRZV,'omitnan');

datStruc.ulckDI20 = mean(datStruc.lckDI(end-19:end),'omitnan');
datStruc.uPreRZV20 = mean(datStruc.preRZV(end-19:end,:),'omitnan');

if saveFlag
    fsave(lickDIFig,[root.name '_lickDI'])
    saveas(vXCorFig,[root.name '_velXCorr'], 'png')
    save([root.name '_datStruc'],'datStruc')
end

%% Get real unit and epoch parameters

datStruc.trueSI = zeros(size(root.good));
datStruc.truePk = zeros(size(root.good));
datStruc.trueLc = zeros(size(root.good));

for i = 1:nUnits
    cc = root.good(i);
    lfpInd = root.info.shankID(root.info.cluster_id == cc)+1;   % Account for 0-indexing
    [datStruc.trueSI(i),~,datStruc.truePk(i),datStruc.trueLc(i),~,~,datStruc.posfr(i,:),datStruc.binedges] = get_SI(root,cc,sess,dbnsz);
    datStruc.frStandRun(i,:) = get_frStandVRun(root,cc,sess);
    [~,~,datStruc.trueVelMdl(i)] = plot_frXvel(root,cc,sess,2,0);
    [datStruc.thetastats(i),datStruc.thetafr(i,:)] = plot_thetaMod(root,cc,lfpInd,2*pi/36,0);
    datStruc.swrfr(i,:) = plot_frXripple(root,cc,sess,ripRef,wlen,histoBnsz,0);
    [~,datStruc.rwdfr(i,:),datStruc.trueRI(i)] = plot_frXrwdtime(root,cc,sess,0.25,5,0);
    [~,datStruc.frMap(:,:,i),datStruc.spkMap(:,:,i)] = get_frXpos(root,cc,sess,0.05,1.85,1);
end

datStruc = subEpochSI(sess,root,datStruc,10);

datStruc.ripRate = size(root.ripStruc(ripRef).ripples,1) / (sum(not(sess.runInds)) / sess.samprate);    %Normalize based on standing periods
datStruc.binpos = datStruc.binedges(1:end-1)+0.5*dbnsz;
datStruc.rwdBin = find(datStruc.binpos > r1pos,1);

%% Find place cells and SPW-R modulation pre and post shift

tic
nShufs = 250;

datStruc.shufSI = zeros(nUnits,nShufs);
datStruc.shufSPWR = zeros(nShufs,nUnits,length(-wlen:histoBnsz:wlen)-1);

for j = 1:nShufs
    [shiftRoot,shiftSess,~] = shiftTrain(root,sess);

    if mod(j,50) == 0
        disp(['Shuffle # ' num2str(j)])
        toc 
    end

    for i = 1:nUnits
        cc = root.good(i);

        [datStruc.shufSI(i,j)] = get_SI(shiftRoot,cc,shiftSess);
        datStruc.shufSPWR(j,i,:) = plot_frXripple(shiftRoot,cc,shiftSess,ripRef,wlen,histoBnsz,0);
        % [~,~,datStruc.shufRI(i,j)] = plot_frXrwdtime(shiftRoot,cc,shiftSess,0.25,5,0);
    end
end

toc

datStruc.sig = sum(datStruc.shufSI >= datStruc.trueSI,2) / nShufs;

disp(['SI p <= 0.05 ' num2str(sum(datStruc.sig <= 0.05)) ' of ' num2str(nUnits) ' units'])
 
%% SPWRs pre and post

for i = 1:nUnits
    cc = root.good(i);
    tmpRipParticip = plot_frXripple(root,cc,sess,ripRef,wlen,histoBnsz,0);
    datStruc.ripParticipation(i) = get_RipParticipation(root,cc,sess,ripRef,wlen);
    p_conf = get_confband(squeeze(datStruc.shufSPWR(:,i,:)),tmpRipParticip);
    datStruc.ripModBinCt(i) = sum(p_conf);
end

disp(['SPWR Modulated p <= 0.05 ' num2str(sum(datStruc.ripModBinCt > 1)) ' of ' num2str(nUnits) ' units'])

%% Get relevant units and waterfall
lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = lyrUnits & hiFRUnits & root.info.uType(root.goodind);
datStruc.useUnits = useUnits;
siUnits = useUnits & datStruc.sig <= 0.05;

r1posInd = find(binpos > r1pos,1);
r2posInd = find(binpos > r2pos,1);
nBins = length(binpos);

[si_wfl_fig,si_wfl_map] = plot_unitsXpos(root,sess,root.good(useUnits));
plot([datStruc.rwdBin datStruc.rwdBin],[0 length(useUnits)+1],'r--','LineWidth',2)
title('Familiar RZ')

% Circularly shifted version
[si_wfl_shft_fig,tmpMap,sortPre] = plot_unitWaterfall(circshift(si_wfl_map,round(nBins/2)-1,2),datStruc.binedges,0,1,0);
xticks([1, round(nBins)/2, nBins]); xticklabels([-90, 0, 90])
plot([datStruc.rwdBin datStruc.rwdBin]+round(nBins/2)-1,[0 length(useUnits)+1],'r--','LineWidth',2)
% plot([r2posInd r2posInd],[0 sum(siBothID)],'k--','LineWidth',2); 
title('Familiar RZ, Centered'); xlabel('Track Position (cm)')

[si_wfl_shft_histo,tmpPks] = plot_unitPkHisto(tmpMap,datStruc.binedges*100,1);
xticks([datStruc.binedges(1),datStruc.binedges(round(nBins/2)),datStruc.binedges(end)]*100); xticklabels([-90, 0, 90])
plot([r2pos r2pos]*100,[0 0.11],'k--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.1])
[ps.lc_pkUniformity, stats.lc_pkUniformity] = pkChi2(tmpPks,datStruc.binedges);
text2bar(si_wfl_shft_histo,'',ps.lc_pkUniformity)

if saveFlag
    fsave(si_wfl_fig,[root.name '_SI_wfl'])
    fsave(si_wfl_shft_fig,[root.name '_SI_shiftwfl'])
    fsave(si_wfl_shft_histo,[root.name '_SI_shiftwfl_histo'])
end

%% PV analysis through time
nBins = length(binpos);

% Compare identity to off-diagonal
idMat = logical(eye(nBins));
% dgMat = logical(spdiags([1 1],[-round(nBins/2) round(nBins/2)],nBins,nBins));

datStruc.pvXlap = get_pvXtime(datStruc.posfr,datStruc.frMap,siUnits,idMat);

pvXtimeFig = figure; hold on;
plot(1:size(datStruc.frMap,1),datStruc.pvXlap)
plot(1:size(datStruc.frMap,1),smoothdata(datStruc.pvXlap),'b','LineWidth',2)
ylim([0 1]); xlabel('Lap'); ylabel('Mean PV Corr. on-diagonal')
legend('Fam->Fam','Smooth')
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

if saveFlag
    fsave(pvXtimeFig,[root.name '_pvXlap'])
end

%% Save for later

if saveFlag
    save([root.name '_datStruc'],'datStruc')
end

%% Prepare new group data
parentDir = "D:\Data\Kelton\analyses\group_analyses"; 
datT = import_xldat(parentDir,"dat_include.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\restCompare';
cd(groupSDir)

cleanInds = datT.include ~= 1;
datT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;
sbase = 'subRestComp_';
fname = 'subRestComp_data1';

dbnsz = 0.05;
histoBnsz = 5;
binedges = 0:5:185;
binpos = 0.025:dbnsz:1.825;
wlen = 150;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

clear ps stats

%% Combine data

region = 'sub';

frDat = [];     % [standFR, runFR]
recID = [];     % [mouseID, recDay, recUnit ID, dist2center, dist2border]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
lcDat = [];     % [si_p, si, pkFR, pkLoc]
lcMap = [];     % [posfr];
siStr = [];     % Struc containing binned SI data per 10 laps
vlDat = [];     % [sig., slope, R]
thDat = [];     % [p, mrl, ang]
thMap = [];     % [thetafr];
rpDat = [];     % [ripParticip, ripModBin]
rpRat = [];     % [ripRate];
rpMap = [];     % [swrfr];
bvDat = [];     % [lckDI, preRZV];
pvStr = [];     % Struc containing binned pv data per laps

ct    = 1;

for i = 1:height(datT)
    if datT.wait(i) ~= 20
        continue
    end

    % === Load data ===
    cd(datT.fpath{i})

    rootfile = dir("*_root.mat");
    try
        load(rootfile.name)
    catch
        disp(['No root file for ' datT.fpath{i}])
    end
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    if datT.sess_type(i) == 2 || datT.sess_type(i) == 3
        datfile = dir("*_RwdShift_Data2.mat");
    elseif datT.sess_type(i) == 1
        datfile = dir("*_datStruc.mat");
    end
    try
        load(datfile.name)
    catch
        disp(['No data file found for ' root.name])
        continue
    end
    disp(root.name)

    if datT.sess_type(i) == 2 || datT.sess_type(i) == 3
        datStruc = frstHalf;
    end

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
        disp(['No ROI units for ' root.name])
        continue
    end

    pkExclude = datStruc.truePk < 0.1;
    useUnits = root.info.lyrID(root.goodind) == 1 & root.info.fr(root.goodind) > 0.1 ...
        & root.info.uType(root.goodind) & roiUnits(root.goodind) & ~pkExclude;
    nCCs = length(root.good);

    % === Concatenate recording data ===
    datStruc.d2cs = get_dist2lyrcenter(root);
    datStruc.d2cs = datStruc.d2cs(root.goodind);
    borderDists = [datT.sh1dist(i) datT.sh2dist(i) datT.sh3dist(i) datT.sh4dist(i)];
    root = get_dist2border(root,borderDists);
    datStruc.d2bs = root.info.bdrDist(root.goodind);
    useCC = logical([useCC; useUnits]);
    recID = [recID; str2num(datT.mouse{i}(end-2:end))*ones(nCCs,1), ...
        datT.session(i)*ones(nCCs,1), root.good, datStruc.d2cs, datStruc.d2bs];

    % === Concatenate Velocity-FR data ===
    for j = 1:length(root.good)
        vlDat = [vlDat; datStruc.trueVelMdl(j).p, datStruc.trueVelMdl(j).b, datStruc.trueVelMdl(j).r];
    end

    % === Concatenate FR data ===
    frDat = [frDat; datStruc.frStandRun];

    % === Concatenate Theta Modulation data ===
    [thAng1, thMRL1, thP1] = get_thAng(datStruc.thetastats);
    thDat = [thDat; thP1', thMRL1', thAng1'];
    thMap = [thMap; datStruc.thetafr];

    % === Concatenate SI and Peak data ===
    lcDat = [lcDat; datStruc.sig, datStruc.trueSI, datStruc.truePk, datStruc.trueLc];
    lcMap = [lcMap; datStruc.posfr];

    % === Concatenate SPWR Modulation data ===
    rpDat = [rpDat; datStruc.ripParticipation', datStruc.ripModBinCt'];
    rpRat = [rpRat; datStruc.ripRate];
    rpMap = [rpMap; datStruc.swrfr];

    % === Concatenate SI over 10-trial blocks ===
    siStr(ct).blockSI = datStruc.subEpochSI;

    % === Concatenate PV over laps ===
    pvStr(ct).blockPV = datStruc.pvXlap;

    % === Concatenate Behavior data ===
    bvDat(ct).lckDI = datStruc.lckDI;
    bvDat(ct).uLckDI = datStruc.ulckDI;
    bvDat(ct).rzVel = datStruc.preRZV;
    bvDat(ct).uRZVel = datStruc.uPreRZV;
    bvDat(ct).uRZVel20 = datStruc.uPreRZV20;
    bvDat(ct).uLckDI20 = datStruc.ulckDI20;
    % bvDat(ct).vMap = frstHalf.spVelMap;

    ct = ct + 1;
end

cd(groupSDir)

save(fname,'frDat','recID','useCC','siStr','lcDat','lcMap','thDat','thMap','rpDat','rpRat','rpMap','vlDat','bvDat','pvStr')

%% Compare old data to new

% Sub Comparisons
% load('D:\Data\Kelton\analyses\ZM002\ZM002_10012025_rec_D5_LLat2\ZM002_10012025_rec_D5_LLat2_RwdShift_Data2.mat');
load('D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\bigcohort\subRwdShift_data5.mat');
r10 = struct('frDat',frDat,'recID',recID,'useCC',useCC,'siStr',siStr,'lcDat',lcDat,'lcMap',lcMap,...
    'thDat',thDat,'thMap',thMap,'rpDat',rpDat,'rpRat',rpRat,'rpMap',rpMap,'vlDat',vlDat,'bvDat',bvDat,'pvStr',pvStr);
load('D:\Data\Kelton\analyses\group_analyses\restCompare\subRestComp_data1.mat');
r20 = struct('frDat',frDat,'recID',recID,'useCC',useCC,'siStr',siStr,'lcDat',lcDat,'lcMap',lcMap,...
    'thDat',thDat,'thMap',thMap,'rpDat',rpDat,'rpRat',rpRat,'rpMap',rpMap,'vlDat',vlDat,'bvDat',bvDat,'pvStr',pvStr);

% CA1 comparisons
% load('D:\Data\Kelton\analyses\ZM001\ZM001_10012025_rec_D5_LLat2\ZM001_10012025_rec_D5_LLat2_RwdShift_Data2.mat');
% load('D:\Data\Kelton\analyses\ZM002\ZM002_09262025_rec_D2_RLat2\ZM002_09262025_rec_D2_RLat2_RwdShift_Data2.mat');
% load('D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\ca1shift\ca1RwdShift_data.mat');

sbase = 'subCompare_';
% cd('D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\bigcohort\compareRest')
% cd('D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\ca1shift\compareRest_ZM001')
% cd('D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\ca1shift\compareRest_ZM002')

ps = struct;
stats = struct;

[~,ps.FRstand,~,stats.FRstand] = ttest2(r10.frDat(r10.useCC,1), r20.frDat(r20.useCC,1));
frStandFig = plotBar2(r10.frDat(r10.useCC,1), r20.frDat(r20.useCC,1));
xticklabels({'10m Rest', '20m Rest'}); text2bar(frStandFig,'FR Stand (Hz)',ps.FRstand);

[~,ps.FRrun,~,stats.FRrun] = ttest2(r10.frDat(r10.useCC,2), r20.frDat(r20.useCC,2));
frRunFig = plotBar2(r10.frDat(r10.useCC,2), r20.frDat(r20.useCC,2));
xticklabels({'10m Rest', '20m Rest'}); text2bar(frRunFig,'FR Run (Hz)',ps.FRrun);

[~,ps.si,~,stats.si] = ttest2(r10.lcDat(r10.useCC & r10.lcDat(:,1) < 0.05,2), r20.lcDat(r20.useCC & r20.lcDat(:,1) < 0.05,2));
siFig = plotBar2(r10.lcDat(r10.useCC & r10.lcDat(:,1) < 0.05,2), r20.lcDat(r20.useCC & r20.lcDat(:,1) < 0.05,2));
xticklabels({'10m Rest', '20m Rest'}); text2bar(siFig,'Spatial Info.',ps.si);

[~,ps.pkFR,~,stats.pkFR] = ttest2(r10.lcDat(r10.useCC & r10.lcDat(:,1) < 0.05,3), r20.lcDat(r20.useCC & r20.lcDat(:,1) < 0.05,3));
pkFRFig = plotBar2(r10.lcDat(r10.useCC & r10.lcDat(:,1) < 0.05,3), r20.lcDat(r20.useCC & r20.lcDat(:,1) < 0.05,3));
xticklabels({'10m Rest', '20m Rest'}); text2bar(pkFRFig,'Peak FR (Hz)',ps.pkFR);

% [thAng, thMRL, thP] = get_thAng(frstHalf.thetastats);
[~,ps.thMRL,~,stats.thMRL] = ttest2(r10.thDat(r10.useCC & r10.thDat(:,1) < 0.05,2), r20.thDat(r20.useCC & r20.thDat(:,1) < 0.05,2));
thMRLFig = plotBar2(r10.thDat(r10.useCC & r10.thDat(:,1) < 0.05,2), r20.thDat(r20.useCC & r20.thDat(:,1) < 0.05,2));
xticklabels({'10m Rest', '20m Rest'}); text2bar(thMRLFig,'Theta MRL',ps.thMRL);

[~,ps.thAng,~,stats.thAng] = ttest2(r10.thDat(r10.useCC & r10.thDat(:,1) < 0.05,3), r20.thDat(r20.useCC & r20.thDat(:,1) < 0.05,3));
thAngFig = plotBar2(r10.thDat(r10.useCC & r10.thDat(:,1) < 0.05,3), r20.thDat(r20.useCC & r20.thDat(:,1) < 0.05,3));
xticklabels({'10m Rest', '20m Rest'}); ylim([-180 180]); text2bar(thAngFig,'Theta Angle',ps.thAng);

% vlSesDat = [];
% for j = 1:length(root.good)
%     vlSesDat = [vlSesDat; frstHalf.trueVelMdl(j).p, frstHalf.trueVelMdl(j).b, frstHalf.trueVelMdl(j).r];
% end
[~,ps.vlb,~,stats.vlb] = ttest2(r10.vlDat(r10.useCC & r10.vlDat(:,1) < 0.05,2), r20.vlDat(r20.useCC & r20.vlDat(:,1) < 0.05,2));
vlbFig = plotBar2(r10.vlDat(r10.useCC & r10.vlDat(:,1) < 0.05,2), r20.vlDat(r20.useCC & r20.vlDat(:,1) < 0.05,2));
xticklabels({'10m Rest', '20m Rest'}); ylim([-1 1]); text2bar(vlbFig,'Velocity Slope',ps.vlb);

[~,ps.vlr,~,stats.vlr] = ttest2(r10.vlDat(r10.useCC & r10.vlDat(:,1) < 0.05,3), r20.vlDat(r20.useCC & r20.vlDat(:,1) < 0.05,3));
vlrFig = plotBar2(r10.vlDat(r10.useCC & r10.vlDat(:,1) < 0.05,3), r20.vlDat(r20.useCC & r20.vlDat(:,1) < 0.05,3));
xticklabels({'10m Rest', '20m Rest'}); text2bar(vlrFig,'Velocity R^2',ps.vlr);

[~,ps.swrPart,~,stats.swrPart] = ttest2(r10.rpDat(r10.useCC & r10.rpDat(:,2) > 1,1), r20.rpDat(r20.useCC & r20.rpDat(:,2) > 1,1));
swrPFig = plotBar2(r10.rpDat(r10.useCC & r10.rpDat(:,2) > 1,1), r20.rpDat(r20.useCC & r20.rpDat(:,2) > 1,1));
xticklabels({'10m Rest', '20m Rest'}); text2bar(swrPFig,'SPWR Particip. Rate',ps.swrPart);

if saveFlag
    saveas(frStandFig,[sbase 'frStand'],'png')
    saveas(frRunFig,[sbase 'frRun'],'png')
    saveas(siFig,[sbase 'spatialinfo'],'png')
    saveas(pkFRFig,[sbase 'peakFieldFR'],'png')
    saveas(thMRLFig,[sbase 'thetaMRL'],'png')
    saveas(thAngFig,[sbase 'thetaAngle'],'png')
    saveas(vlbFig,[sbase 'velocityB'],'png')
    saveas(vlrFig,[sbase 'velocityR'],'png')
    saveas(swrPFig,[sbase 'spwrParticipation'],'png')
end
