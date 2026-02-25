% EC
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05062025_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05072025_rec_D2_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05082025_rec_D3_RLat2';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05092025_rec_D4_RMed2';
% spath = 'D:\Data\Kelton\analyses\KW054\KW054_07292025_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW054\KW054_07302025_rec_D2_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW054\KW054_07312025_rec_D3_RMed2';
% spath = 'D:\Data\Kelton\analyses\KW084\KW084_12172025_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW084\KW084_12182025_rec_D2_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW083\KW083_12172025_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW083\KW083_12182025_rec_D2_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW083\KW083_12192025_rec_D3_RLat2';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_ec_pilot_dat.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

epochInd = round(sess.nlaps/2);
nUnits = length(root.good);
saveFlag = 1;

dbnsz = 0.05;
histoBnsz = 5;
wlen = 150;
ripRef = root.ripRef;
r1pos = 0.9;    % 90 cm
binpos = 0.025:dbnsz:1.825;
tmpD = root.info.depth(root.goodind);

sbase = root.name; 

%% Plot Lick DI behavior

lickDIFig = plot_lickDiscrim(sess,[r1pos 1.8]*100,30,10);

if saveFlag
    saveas(lickDIFig,[sbase '_lickDI'], 'png')
end

%% Get real unit and epoch parameters

clear datStruc
datStruc.name   = root.name;
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
    % datStruc.swrfr(i,:) = plot_frXripple(root,cc,sess,ripRef,wlen,histoBnsz,0);
    [~,datStruc.rwdfr(i,:),datStruc.trueRI(i)] = plot_frXrwdtime(root,cc,sess,0.25,5,0);
    [datStruc.optobins,datStruc.optoMat(i,:,:)] = plot_frXopto(root,cc,sess,0.002,0.5,0);
    tmp = get_firstPk(datStruc,i);
    datStruc.optoPkT(i) = tmp;
end

datStruc = subEpochSI(sess,root,datStruc,10);

% datStruc.ripRate = size(root.ripStruc(ripRef).ripples,1) / (sum(not(sess.runInds)) / sess.samprate);    %Normalize based on standing periods
datStruc.binpos = datStruc.binedges(1:end-1)+0.5*dbnsz;

%% Find place cells and SPW-R modulation pre and post shift
% Methods: 
% Kitanishi et al. 2021: Exceed 99th percentile of SI shuffle
% Grienberger & Magee 2022:
%   Field: contiguous area with >= 20% peak FR
%   Induction: First lap with >3SD of in-field FR, above out-field noise
%              and activity in 2/5 subsequent laps also >3SD
%   Stability: Significant >3SD activity in 30% of post-induction laps

% Add method for jittering opto tagging

tic
nShufs = 250;

datStruc.shufSI = zeros(nUnits,nShufs);
datStruc.shufSPWR = zeros(nShufs,nUnits,length(-wlen:histoBnsz:wlen)-1);

for j = 1:nShufs
    [shiftroot,shiftsess,shiftInd] = shiftTrain(root,sess);

    if mod(j,50) == 0
        disp(['Shuffle # ' num2str(j)])
        toc 
    end

    for i = 1:nUnits
        cc = root.good(i);

        [datStruc.shufSI(i,j)] = get_SI(shiftroot,cc,shiftsess);

        % frstHalf.shufSPWR(j,i,:) = plot_frXripple(shiftFrst,cc,shiftSessFrst,ripRef,wlen,histoBnsz,0);

        % [~,~,frstHalf.shufRI(i,j)] = plot_frXrwdtime(shiftFrst,cc,shiftSessFrst,0.25,5,0);
        % [~,~,lastHalf.shufRI(i,j)] = plot_frXrwdtime(shiftLast,cc,shiftSessLast,0.25,5,0);
    end
end

toc

datStruc.sig = sum(datStruc.shufSI > datStruc.trueSI,2) / nShufs;

disp(['Sig SI p <= 0.05 ' num2str(sum(datStruc.sig <= 0.05)) ' of ' num2str(nUnits) ' units'])

%% Save for later

if saveFlag
    save([root.name '_ec_pilot_dat'],'datStruc')
end

%% Plot unit parameters by depth

for i = 1:length(root.good)
    tmp = get_firstPk(datStruc,i);
    datStruc.optoPkT(i) = tmp;
end

thDepthFig = plot_datXdepth(root,datStruc,1,0,0,2); % theta
siDepthFig = plot_datXdepth(root,datStruc,1,0,0,3); % spatial info
opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time

if saveFlag
    saveas(thDepthFig,[sbase '_ec_pilot_thXdepth.png'])
    saveas(siDepthFig,[sbase '_ec_pilot_siXdepth.png'])
    saveas(opDepthFig,[sbase '_ec_pilot_opXdepth.png'])
end

%% Splitting by gamma and depth
% Fernandez-Ruiz et al., 2021 Gamma-S = 30-50; Gamma-M = 60-80; Gamma-L = 100-150

root = addRootLFP(spath);
root = get_lfpXdepth(root,[6 10; 150 250; 30 50; 60 80; 100 150]);
tmpLFPDensityFig = plot_lfpXdepth(root);
ylim([0 max(root.lfpinfo.lfpDepth)])
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])

saveas(tmpLFPDensityFig,[sbase '_LFPPower.png'])

%% Assign putative layer

for i = 1:length(root.good)
    datStruc.thAng(i) = rad2deg(datStruc.thetastats(i).ang);
    datStruc.thSig(i) = datStruc.thetastats(i).p < 0.05;
end

shortLat = datStruc.optoPkT <= 0.01;

% figure; hold on
% plot(datStruc.thAng+180,datStruc.optoPkT,'ko')
% plot(datStruc.thAng(shortLat)+180,datStruc.optoPkT(shortLat),'r*')

% PCA attempt
% clusterDat = [normalize(datStruc.thAng'+180,'range'),normalize(datStruc.optoPkT','range'),normalize(tmpD,'range')];
% 
% [coeff,score,latent,tsquare,explained] = pca(clusterDat);
% 
% figure; 
% plot3(score(:,1),score(:,2),score(:,3),'k.')

% moving average of opto pulse plus low latency units
[sortDepth,sortInd] = sort(tmpD);
sortOpto = datStruc.optoPkT(sortInd);
meanOpPk = smoothdata(sortOpto,'movmean',10);
ec3stt = sortDepth(find(meanOpPk < 0.025,1,'first'));
ec3end = sortDepth(find(meanOpPk < 0.025,1,'last'));

figure; hold on
plot(sortOpto)
plot(meanOpPk)
xlabel('unit (sorted by depth)')
ylabel('First opto peak post-pulse (s)')
legend('Raw','Moving Mean')
set(gca,'FontSize',12,'FontName','Arial')

ec5inds = tmpD > ec3end;
ec2inds = tmpD < ec3stt;
ec3inds = (tmpD >= ec3stt & tmpD <= ec3end) | datStruc.optoPkT' < 0.015;

datStruc.lyrID = zeros(size(datStruc.trueSI));
datStruc.lyrID(ec5inds) = 5;
datStruc.lyrID(ec2inds) = 2;
datStruc.lyrID(ec3inds) = 3;

opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time
plot(datStruc.optoPkT(datStruc.lyrID == 5)*1000, tmpD(datStruc.lyrID == 5),'k*')
plot(datStruc.optoPkT(datStruc.lyrID == 3)*1000, tmpD(datStruc.lyrID == 3),'c*')
plot(datStruc.optoPkT(datStruc.lyrID == 2)*1000, tmpD(datStruc.lyrID == 2),'m*')
% legend('Data','EC5','EC3','EC2')
set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]); xlim([0 100]); xticks([0 100]); xticklabels([0 100])

thDepthFig = plot_datXdepth(root,datStruc,1,0,0,2); % theta
plot(datStruc.thAng(datStruc.lyrID == 5)+180, tmpD(datStruc.lyrID == 5),'k*')
plot(datStruc.thAng(datStruc.lyrID == 3)+180, tmpD(datStruc.lyrID == 3),'c*')
plot(datStruc.thAng(datStruc.lyrID == 2)+180, tmpD(datStruc.lyrID == 2),'m*')
% legend('Data','EC5','EC3','EC2')
set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]);

saveas(opDepthFig,[sbase '_opXdepth.png'])
saveas(thDepthFig,[sbase '_thXdepth.png'])

%% Plot theta filtered lfp by depth over time window
thlfp = bandpass(root.lfp', [6 10], root.fs_lfp);
twin = [600 630];
indwin = twin*root.fs_lfp ;

thXdepthFig = figure; hold on
set(gcf,'units','normalized','position',[0.45 0.05 0.14 0.6])
for i = 1:size(thlfp,2)
    plot(sess.ts(root.lfp_tsb(indwin(1):indwin(2))), 100*thlfp(indwin(1):indwin(2),i) + root.lfpinfo.lfpDepth(i),'k')
end
ylim([-30 2600]); xlim([615 616]); xlabel('Time (s)'); ylabel('Depth (um)')
set(gca,'FontSize',12,'FontName','Arial')
grid on
saveas(thXdepthFig,[sbase '_thLFPXdepth.png'])

%% Descriptive statistics and graphs
% root.good(sigFrst <= 0.05 & sigLast <= 0.05 & root.info.fr(root.goodind)' > 0.1 & root.info.lyrID(root.goodind)' == 1)

% hiFRUnits = root.info.fr(root.goodind) > 0.1 & datStruc.truePk > 0.1;
% useUnits = hiFRUnits & root.info.uType(root.goodind);   % hi FR and not IN
% siUnits = useUnits & datStruc.sig < 0.05;
% 
% nUseUnits = sum(useUnits);
% 
% xcoords = ones(nUseUnits,1);
% shiftbins = [-fliplr(datStruc.binedges) datStruc.binedges(2:end)];
% 
% % Firing rate stand vs run
% [~,ps.standFR_firstlast] = ttest(datStruc.frStandRun(useUnits,1),lastHalf.frStandRun(useUnits,1));
% [~,ps.runFR_firstlast] = ttest(datStruc.frStandRun(useUnits,2),lastHalf.frStandRun(useUnits,2));
% 
% % Spatial Information
% [~,ps.SI_firstlast] = ttest(datStruc.trueSI(siUnits),lastHalf.trueSI(siUnits));
% 
% lcBothSIFig = plotBar2(datStruc.trueSI(siUnits), lastHalf.trueSI(siUnits));
% ylabel('Spatial Information (Bits/spike)')
% 
% if saveFlag
%     sbase = root.name;
%     saveas(frStndRunFig,[sbase '_RwdShift_FRStandRun.png'])
%     saveas(lcBothSIFig,[sbase '_RwdShift_SI.png'])
% end

%% Peak distribution - all EC
datStruc.rwdBin = find(datStruc.binpos > r1pos,1);

ecFieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(siUnits)),datStruc.binedges);

[spWflFig] = plot_unitsXpos(root,sess,root.good(siUnits));
plot([datStruc.rwdBin datStruc.rwdBin],[0 length(siUnits)+1],'r--','LineWidth',2)
% plot(normalize(ecFieldDistro,'range').*sum(siUnits)./5,'k','LineWidth',1)
title(replace(root.name,"_"," "))

peakDistroFig = plotDistroHisto(ecFieldDistro,binpos,r1pos);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')

if saveFlag
    saveas(spWflFig,[sbase '_ec_pilot_waterfall.png'])
    saveas(peakDistroFig,[sbase '_ec_pilot_PeakDistro.png'])
end

%% Peak distribution - by EC layer
nBins = length(datStruc.binpos);

[wflFig2,tmpMap2,sortEC2] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 2 & siUnits,:),datStruc.binedges);
ec2FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 2 & siUnits)),datStruc.binedges);
% peak2DistroFig = plotDistroHisto(ec2FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

[wflFig3,tmpMap3,sortEC3] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 3 & siUnits,:),datStruc.binedges);
ec3FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 3 & siUnits)),datStruc.binedges);
% peak3DistroFig = plotDistroHisto(ec3FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

[wflFig5,tmpMap5,sortEC5] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 5 & siUnits,:),datStruc.binedges);
ec5FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 5 & siUnits)),datStruc.binedges);
% peak5DistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

% Stack all together
ecSortMap = [tmpMap2; tmpMap3; tmpMap5];
nUnits = [size(tmpMap5,1) size(tmpMap3,1) size(tmpMap2,1)];
wflFigLyrs = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.05 0.25 0.8])
imagesc(ecSortMap,[prctile(ecSortMap,1,'all'), prctile(ecSortMap,98,'all')]);
plot([0.5 nBins],[nUnits(3) nUnits(3)],'k--','LineWidth',2)
plot([0.5 nBins],[sum(nUnits(2:3)) sum(nUnits(2:3))],'k--','LineWidth',2)
plot([datStruc.rwdBin datStruc.rwdBin],[0 sum(nUnits)],'w--','LineWidth',2)
colormap("parula")
cbar = colorbar; clim([0 0.98]);
xticks(1:10:nBins+1)
xticklabels(100*datStruc.binedges(1:10:end))
xlim([0.15 nBins+0.85]); ylim([0.15 sum(nUnits)+0.85])
ylabel('Unit #'); ylabel(cbar,'FR (Normalized)','FontSize',12,'Rotation',90)
xlabel('Track Position (cm)')
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

% Plot distros together
peakLyrDistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
plot(binpos*100,ec2FieldDistro./sum(ec2FieldDistro),'Color','m','LineWidth',2);
plot(binpos*100,ec3FieldDistro./sum(ec3FieldDistro),'Color','c','LineWidth',2);
plot(binpos*100,ec5FieldDistro./sum(ec5FieldDistro),'Color','k','LineWidth',2);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')

if saveFlag
    % saveas(wflFig2,[sbase '_ec_pilot_waterfall_EC2.png'])
    % saveas(wflFig3,[sbase '_ec_pilot_waterfall_EC3.png'])
    % saveas(wflFig5,[sbase '_ec_pilot_waterfall_EC5.png'])
    saveas(wflFigLyrs,[sbase '_ec_pilot_waterfall_EC2-5.png'])
    saveas(peakLyrDistroFig,[sbase '_ec_pilot_waterfall_EC2-5_histo.png'])
end

%% Putative 1D grid cells
% From Wen et al., 2024
% Z-score spatial FR for each unit
% Compute spectrogram using 16lap / 30m window
% Look for 3 peak structure in PSD

for i = 1:length(root.good)
    cc = root.good(i);

    dbnsz = 0.03;
    [tmpedges,tmpbinfr] = get_frXpos(root,cc,sess,dbnsz);

    spFR = reshape(tmpbinfr',[numel(tmpbinfr),1]);
    spFR = normalize(spFR);

    wn = (length(tmpedges)-1) * 16;     % 16 laps worth of segments
    ovlp = round(wn*0.875);
    [pxx,f] = pwelch(spFR,wn,ovlp,10000,3);

    figure; plot(f,pxx)
    ylim([0 50])
    ylabel('PSD'); xlabel('Freq')
end

%% Epoch into pre/post reward shift for reward shift sessions

epoch1Ind = 30;
frstHalfInds   = [sess.ind(1) sess.lapend(epoch1Ind)];
lastHalfInds   = [sess.lapstt(end-epoch1Ind) sess.ind(end)];
[sess1, root1] = epochStruc(sess,root,frstHalfInds);
sess1.valTrials = 2:sess.nlaps;
[sess2, root2] = epochStruc(sess,root,lastHalfInds);

% Get real unit and epoch parameters

clear epoch1 epoch2
epoch1.name   = root.name;
epoch2.name   = root.name;

for i = 1:length(root.good)
    cc = root.good(i);
    [~,~,~,epoch1.trueLc(i),~,~,epoch1.posfr(i,:)] = get_SI(root1,cc,sess1,dbnsz);
    [~,~,~,epoch2.trueLc(i),~,~,epoch2.posfr(i,:)] = get_SI(root2,cc,sess2,dbnsz);
end

% epoch1 = subEpochSI(sess1,root1,epoch1,10);
% epoch2 = subEpochSI(sess2,root2,epoch2,10);

%% Plot epoch waterfalls
[wflFig1_EC2,tmpMap1_EC2] = plot_unitWaterfall(epoch1.posfr(datStruc.lyrID == 2 & siUnits,:),datStruc.binedges,sortEC2);
ec2FieldDistro1 = histcounts(datStruc.binpos(epoch1.trueLc(datStruc.lyrID == 2 & siUnits)),datStruc.binedges);
[wflFig2_EC2,tmpMap2_EC2] = plot_unitWaterfall(epoch2.posfr(datStruc.lyrID == 2 & siUnits,:),datStruc.binedges,sortEC2);
ec2FieldDistro2 = histcounts(datStruc.binpos(epoch2.trueLc(datStruc.lyrID == 2 & siUnits)),datStruc.binedges);
% Plot distros together
peakEC2DistroFig = plotDistroHisto(ec2FieldDistro1,binpos,r1pos);
plot(binpos*100,ec2FieldDistro2./sum(ec2FieldDistro2),'Color','m','LineWidth',2);
legend({'First half','Last half'},'Location','northwest')
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])

[wflFig1_EC3,tmpMap1_EC3] = plot_unitWaterfall(epoch1.posfr(datStruc.lyrID == 3 & siUnits,:),datStruc.binedges,sortEC3);
ec3FieldDistro1 = histcounts(datStruc.binpos(epoch1.trueLc(datStruc.lyrID == 3 & siUnits)),datStruc.binedges);
[wflFig2_EC3,tmpMap2_EC3] = plot_unitWaterfall(epoch2.posfr(datStruc.lyrID == 3 & siUnits,:),datStruc.binedges,sortEC3);
ec3FieldDistro2 = histcounts(datStruc.binpos(epoch2.trueLc(datStruc.lyrID == 3 & siUnits)),datStruc.binedges);
% Plot distros together
peakEC3DistroFig = plotDistroHisto(ec3FieldDistro1,binpos,r1pos);
plot(binpos*100,ec3FieldDistro2./sum(ec3FieldDistro2),'Color','c','LineWidth',2);
legend({'First half','Last half'},'Location','northwest')
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])

[wflFig1_EC5,tmpMap1_EC5] = plot_unitWaterfall(epoch1.posfr(datStruc.lyrID == 5 & siUnits,:),datStruc.binedges,sortEC5);
ec5FieldDistro1 = histcounts(datStruc.binpos(epoch1.trueLc(datStruc.lyrID == 5 & siUnits)),datStruc.binedges);
[wflFig2_EC5,tmpMap2_EC5] = plot_unitWaterfall(epoch2.posfr(datStruc.lyrID == 5 & siUnits,:),datStruc.binedges,sortEC5);
ec5FieldDistro2 = histcounts(datStruc.binpos(epoch2.trueLc(datStruc.lyrID == 5 & siUnits)),datStruc.binedges);
% Plot distros together
peakEC5DistroFig = plotDistroHisto(ec5FieldDistro1,binpos,r1pos);
plot(binpos*100,ec5FieldDistro2./sum(ec5FieldDistro2),'Color','k','LineWidth',2);
legend({'First half','Last half'},'Location','northwest')
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])

if saveFlag
    saveas(wflFig1_EC2,[sbase '_ec_pilot_waterfall_EC2_e1.svg'])
    saveas(wflFig2_EC2,[sbase '_ec_pilot_waterfall_EC2_e2.svg'])
    saveas(peakEC2DistroFig,[sbase '_ec_pilot_waterfall_EC2_e1_histo.svg'])
    saveas(wflFig1_EC3,[sbase '_ec_pilot_waterfall_EC3_e1.svg'])
    saveas(wflFig2_EC3,[sbase '_ec_pilot_waterfall_EC3_e2.svg'])
    saveas(peakEC3DistroFig,[sbase '_ec_pilot_waterfall_EC3_e1_histo.svg'])
    saveas(wflFig1_EC5,[sbase '_ec_pilot_waterfall_EC5_e1.svg'])
    saveas(wflFig2_EC5,[sbase '_ec_pilot_waterfall_EC5_e2.svg'])
    saveas(peakEC5DistroFig,[sbase '_ec_pilot_waterfall_EC5_e1_histo.svg'])
end

%% Save for later

if saveFlag
    save([root.name '_ec_pilot_dat'],'datStruc','epoch1','epoch2')
end

%% Combine Sessions/Animals

parentDir = "D:\Data\Kelton\analyses\group_analyses"; 
datT = import_xldat(parentDir,"dat_include.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\dp2_progreport';
cd(groupSDir) 

cleanInds = datT.include == 0;
datT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;
sbase = 'dp2_progReport_';

dbnsz = 0.05;
histoBnsz = 5;
binpos = 0.025:dbnsz:1.825;
binedges = 0:0.05:1.85;
wlen = 150;
r1pos = 0.9;    % 10 cm
vColors2 = [0.5 0.5 1; 0.75 0.75 1];
r1posInd = find(binpos > r1pos,1);

%% Load old data
load('D:\Data\Kelton\analyses\group_analyses\dp2_progreport\dp2_progReport_data.mat');

%% 
recID = [];     % [mouseID, recDay, recUnitID, recLayer]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
lcMap = [];     % [frstHalf.posfr, lastHalf.posfr];
siStr = [];     % Struc containing binned SI data per 10 laps
% vlDat = [];     % [frstHalf sig., frstHalf slope, frstHalf R, lastHalf sig., lastHalf slope, lastHalf R]
% thDat = [];     % [frstHalf.p, frstHalf.mrl, frstHalf.ang, lastHalf.p, lastHalf.mrl, lastHalf.ang]
% thMap = [];     % [frstHalf.thetafr, lastHalf.thetafr];
% bvDat = [];     % [frstHalf.lckDI, frstHalf.preRZV, lastHalf.lckDI, lastHalf.preRZV

clear ps stats

ct    = 1;

for i = 1:height(datT)
    if datT.sess_type(i) ~= 6
        continue
    end

    % === Load data ===
    cd(datT.fpath{i})

    rootfile = dir("*_root.mat");
    load(rootfile.name)
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    epochfile = dir("*_ec_pilot_dat.mat");
    load(epochfile.name)
    disp(root.name)

    useUnits = root.info.fr(root.goodind) > 0.1 & root.info.uType(root.goodind);
    nCCs = length(root.good);

    % === Concatenate recording data ===
    useCC = logical([useCC; useUnits]);
    recID = [recID; str2num(datT.mouse{i}(end-2:end))*ones(nCCs,1), datT.session(i)*ones(nCCs,1), root.good, datStruc.lyrID];

    % % === Concatenate Velocity-FR data ===
    % for j = 1:length(root.good)
    %     vlDat = [vlDat; frstHalf.trueVelMdl(j).p, frstHalf.trueVelMdl(j).b, frstHalf.trueVelMdl(j).r...
    %         lastHalf.trueVelMdl(j).p, lastHalf.trueVelMdl(j).b, lastHalf.trueVelMdl(j).r];
    % end
    % 
    % % === Concatenate FR data ===
    % frDat = [frDat; frstHalf.frStandRun, lastHalf.frStandRun];

    % % === Concatenate Theta Modulation data ===
    % for j = 1:length(root.good)
    %     thDat = [thDat; frstHalf.thetastats(j).p, frstHalf.thetastats(j).mrl, frstHalf.thetastats(j).ang,...
    %         lastHalf.thetastats(j).p, lastHalf.thetastats(j).mrl, lastHalf.thetastats(j).ang];
    % end
    % thMap = [thMap; frstHalf.thetafr lastHalf.thetafr];

    % === Concatenate SI and Peak data ===
    lcDat = [lcDat; datStruc.sig, datStruc.trueSI, datStruc.truePk, datStruc.trueLc];
    lcMap = [lcMap; datStruc.posfr];

    % === Concatenate SI over 10-trial blocks ===
    siStr(ct).blockSI = datStruc.subEpochSI;

    % % === Concatenate Behavior data ===
    % bvDat(ct).lckDI = datStruc.lckDI;
    % bvDat(ct).uLckDI = datStruc.ulckDI;
    % bvDat(ct).rzVel = datStruc.preRZV;
    % bvDat(ct).uRZVel = datStruc.uPreRZV;
    
    ct = ct + 1;
end

cd(groupSDir)

save([sbase 'data'],'recID','useCC','lcDat','lcMap','siStr')

%%
nBins = length(binpos);

useEC2 = useCC & recID(:,4) == 2 & lcDat(:,1) < 0.05;
useEC3 = useCC & recID(:,4) == 3 & lcDat(:,1) < 0.05;
useEC5 = useCC & recID(:,4) == 5 & lcDat(:,1) < 0.05;

[wflFig2,tmpMap2,sortEC2] = plot_unitWaterfall(lcMap(useEC2,:),binedges);
ec2FieldDistro = histcounts(binpos(lcDat(useEC2,4)),binedges);
% peak2DistroFig = plotDistroHisto(ec2FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

[wflFig3,tmpMap3,sortEC3] = plot_unitWaterfall(lcMap(useEC3,:),binedges);
ec3FieldDistro = histcounts(binpos(lcDat(useEC3,4)),binedges);
% peak3DistroFig = plotDistroHisto(ec3FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

[wflFig5,tmpMap5,sortEC5] = plot_unitWaterfall(lcMap(useEC5,:),binedges);
ec5FieldDistro = histcounts(binpos(lcDat(useEC5,4)),binedges);
% peak5DistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

% Stack all together
ecSortMap = [tmpMap2; tmpMap3; tmpMap5];
nUnits = [size(tmpMap5,1) size(tmpMap3,1) size(tmpMap2,1)];
wflFigLyrs = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.05 0.25 0.8])
imagesc(ecSortMap,[prctile(ecSortMap,1,'all'), prctile(ecSortMap,98,'all')]);
plot([0.5 nBins],[nUnits(3) nUnits(3)],'k--','LineWidth',2)
plot([0.5 nBins],[sum(nUnits(2:3)) sum(nUnits(2:3))],'k--','LineWidth',2)
plot([r1posInd r1posInd],[0 sum(nUnits)],'w--','LineWidth',2)
colormap("parula")
cbar = colorbar; clim([0 0.98]);
xticks(1:10:nBins+1)
xticklabels(100*binedges(1:10:end))
xlim([0.15 nBins+0.85]); ylim([0.15 sum(nUnits)+0.85])
ylabel('Unit #'); ylabel(cbar,'FR (Normalized)','FontSize',12,'Rotation',90)
xlabel('Track Position (cm)')
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

% Plot distros together
peakLyrDistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
plot(binpos*100,ec2FieldDistro./sum(ec2FieldDistro),'Color','m','LineWidth',2);
plot(binpos*100,ec3FieldDistro./sum(ec3FieldDistro),'Color','c','LineWidth',2);
plot(binpos*100,ec5FieldDistro./sum(ec5FieldDistro),'Color','k','LineWidth',2);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')

if saveFlag
    % saveas(wflFig2,[sbase '_ec_pilot_waterfall_EC2.png'])
    % saveas(wflFig3,[sbase '_ec_pilot_waterfall_EC3.png'])
    % saveas(wflFig5,[sbase '_ec_pilot_waterfall_EC5.png'])
    saveas(wflFigLyrs,[sbase '_ec_pilot_waterfall_EC2-5.svg'])
    saveas(peakLyrDistroFig,[sbase '_ec_pilot_waterfall_EC2-5_histo.svg'])
end

%% Functions


function [fhandle] = plotDistroHisto(distro1,binpos,rzPos)
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

fhandle = figure; hold on
plot(binpos*100,distro1./sum(distro1),'Color',vColors2(1,:),'LineWidth',2);
plot([rzPos rzPos]*100,[0 0.15],'k--','HandleVisibility','off')
ylim([0 max(distro1./sum(distro1),[],'all')])
xlim([binpos(1)*100-1 binpos(end)*100+1])
ylabel('P(field peak)')
set(gca,'FontSize',12,'FontName','Arial')
end
