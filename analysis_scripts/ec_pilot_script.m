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
% spath = 'D:\Data\Kelton\analyses\KW093\KW093_03262026_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW093\KW093_03272026_rec_D2_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW093\KW093_03282026_rec_D3_RLat2';
% spath = 'D:\Data\Kelton\analyses\KW093\KW093_03292026_rec_D4_RMed2';
% spath = 'D:\Data\Kelton\analyses\KW098\KW098_05042026_rec_D1_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW098\KW098_05052026_rec_D2_RLat2';
spath = 'D:\Data\Kelton\analyses\KW098\KW098_05062026_rec_D3_RMed1';

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
    % [~,datStruc.rwdfr(i,:),datStruc.trueRI(i)] = plot_frXrwdtime(root,cc,sess,0.25,5,0);
    [datStruc.optobins,datStruc.optoMat(i,:,:)] = plot_frXopto(root,cc,sess,0.002,0.5,0);
    datStruc.optoPkT(i) = get_firstPk(datStruc,i);
end

% datStruc = subEpochSI(sess,root,datStruc,10);

% datStruc.ripRate = size(root.ripStruc(ripRef).ripples,1) / (sum(not(sess.runInds)) / sess.samprate);    %Normalize based on standing periods
datStruc.binpos = datStruc.binedges(1:end-1)+0.5*dbnsz;

%% Find place cells pre and post shift

nShufs = 250;

datStruc = get_shufParams(root,root.good,sess,datStruc,nShufs);
% tic
% 
% datStruc.shufSI = zeros(nUnits,nShufs);
% % datStruc.shufSPWR = zeros(nShufs,nUnits,length(-wlen:histoBnsz:wlen)-1);
% 
% for j = 1:nShufs
%     [shiftroot,shiftsess,shiftInd] = shiftTrain(root,sess);
% 
%     if mod(j,50) == 0
%         disp(['Shuffle # ' num2str(j)])
%         toc 
%     end
% 
%     for i = 1:nUnits
%         cc = root.good(i);
% 
%         [datStruc.shufSI(i,j)] = get_SI(shiftroot,cc,shiftsess);
% 
%         % frstHalf.shufSPWR(j,i,:) = plot_frXripple(shiftFrst,cc,shiftSessFrst,ripRef,wlen,histoBnsz,0);
% 
%         % [~,~,frstHalf.shufRI(i,j)] = plot_frXrwdtime(shiftFrst,cc,shiftSessFrst,0.25,5,0);
%         % [~,~,lastHalf.shufRI(i,j)] = plot_frXrwdtime(shiftLast,cc,shiftSessLast,0.25,5,0);
%     end
% end
% 
% toc
% 
% datStruc.sigSI = sum(datStruc.shufSI > datStruc.trueSI,2) / nShufs;

disp(['Sig SI p <= 0.05 ' num2str(sum(datStruc.sigSI <= 0.05)) ' of ' num2str(nUnits) ' units'])

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
try
    opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time
catch
end

if saveFlag
    saveas(thDepthFig,[sbase '_ec_pilot_thXdepth.png'])
    saveas(siDepthFig,[sbase '_ec_pilot_siXdepth.png'])
    try
        saveas(opDepthFig,[sbase '_ec_pilot_opXdepth.png'])
    catch
    end
end

%% Splitting by gamma and depth
% Fernandez-Ruiz et al., 2021 Gamma-S = 30-50; Gamma-M = 60-80; Gamma-L = 100-150

root = addRootLFP(spath);
root = get_lfpXdepth(root,[6 10; 150 250; 30 50; 60 80; 100 150]);
tmpLFPDensityFig = plot_lfpXdepth(root);
ylim([0 max(root.lfpinfo.lfpDepth)])
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])

saveas(tmpLFPDensityFig,[sbase '_LFPPower.png'])

%% Assign putative layer - using opto

[datStruc.thAng, ~, thSigs] = get_thAng(datStruc.thetastats);
datStruc.thSig = thSigs < 0.05;

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

%% Assign putative layer - using theta angle

[datStruc.thAng, ~, thSigs] = get_thAng(datStruc.thetastats);
datStruc.thSig = thSigs < 0.05;

ec2inds = datStruc.thAng > -90 & datStruc.thAng < 90;
ec3inds = ~ec2inds; %-180:-90 and 90:180;

dBins = 0:50:max(root.info.depth);
lyrct = histcounts(tmpD,dBins); % All good by depth bin
lyrsg = histcounts(tmpD(datStruc.thSig),dBins); % Significantly theta by depth
lyrtr2 = histcounts(tmpD(ec2inds),dBins); % Theta tuning by depth
ec5stt = dBins(find(smooth(lyrsg./lyrct) <= 0.1,1)); % Estimate L5 start using cutoff of ratio of significant / all
ec3stt = dBins(find(smooth(lyrtr2./lyrct) <= 0.4,1)); % Estimate L3 start using cutoff of ratio of L2 theta / all

ec5inds = tmpD' > ec5stt & ~datStruc.thSig;
ec2inds = ec2inds & tmpD' < ec3stt + 300;
try
    ec3inds = (~ec2inds & ~ec5inds) | datStruc.optoPkT < 0.015;
catch
    ec3inds = ~ec2inds & ~ec5inds;
end

datStruc.lyrID = zeros(size(root.good));
datStruc.lyrID(ec2inds) = 2;
datStruc.lyrID(ec3inds) = 3;
datStruc.lyrID(ec5inds) = 5;

thDepthFig = plot_datXdepth(root,datStruc,1,0,0,2); % theta
plot(datStruc.thAng(datStruc.lyrID == 5)+180, tmpD(datStruc.lyrID == 5),'k*')
plot(datStruc.thAng(datStruc.lyrID == 3)+180, tmpD(datStruc.lyrID == 3),'c*')
plot(datStruc.thAng(datStruc.lyrID == 2)+180, tmpD(datStruc.lyrID == 2),'m*')
xlim([0 360]); ylabel('Distance from tip (mm)'); xlabel('Theta angle')
% legend('Data','EC5','EC3','EC2')
% set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]);
%% Plot minis for RPPR 2026

miniOptoF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.4])
scatter(datStruc.thAng(datStruc.lyrID == 5)+180, tmpD(datStruc.lyrID == 5)/1000,'k', 'filled')
scatter(datStruc.thAng(datStruc.lyrID == 3)+180, tmpD(datStruc.lyrID == 3)/1000,'c', 'filled')
scatter(datStruc.thAng(datStruc.lyrID == 2)+180, tmpD(datStruc.lyrID == 2)/1000,'m', 'filled')
xlim([0 360]); ylim([0 2.5]); ylabel('Distance from tip (mm)'); xlabel('\theta angle')
normThPwr = normalize(root.uPSD(1,:),'range',[0 360]);
plot(normThPwr,root.lfpinfo.lfpDepth/1000,'k')
set(gca,'FontName','Arial','FontSize',12)

miniThetaF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.4])
scatter(datStruc.optoPkT(datStruc.lyrID == 5)*1000, tmpD(datStruc.lyrID == 5)/1000,'k', 'filled')
scatter(datStruc.optoPkT(datStruc.lyrID == 3)*1000, tmpD(datStruc.lyrID == 3)/1000,'c', 'filled')
scatter(datStruc.optoPkT(datStruc.lyrID == 2)*1000, tmpD(datStruc.lyrID == 2)/1000,'m', 'filled')
xlim([0 150]); ylim([0 2.5]); ylabel('Distance from tip (mm)'); xlabel('Post-opto Peak time (ms)')
set(gca,'FontName','Arial','FontSize',12)

if saveFlag
    fsave(miniOptoF,[sbase '_depthXopto_mini'],1,1)
    fsave(miniThetaF,[sbase '_depthXtheta_mini'],1,1)
end

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

hiFRUnits = root.info.fr(root.goodind) > 0.1 & datStruc.truePk > 0.1;
useUnits = hiFRUnits & root.info.uType(root.goodind);   % hi FR and not IN
siUnits = useUnits & datStruc.sigSI < 0.05;

nUseUnits = sum(useUnits);

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

[wflFig2,tmpMap2,sortEC2] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 2 & siUnits,:),datStruc.binedges,0,1,0);
ec2FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 2 & siUnits)),datStruc.binedges);
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.32]); ylabel(''); % xlabel('Track position (cm)')
% peak2DistroFig = plotDistroHisto(ec2FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

[wflFig3,tmpMap3,sortEC3] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 3 & siUnits,:),datStruc.binedges,0,1,0);
ec3FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 3 & siUnits)),datStruc.binedges);
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.32]); ylabel(''); % xlabel('Track position (cm)')
% peak3DistroFig = plotDistroHisto(ec3FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

[wflFig5,tmpMap5,sortEC5] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 5 & siUnits,:),datStruc.binedges,0,1,0);
ec5FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 5 & siUnits)),datStruc.binedges);
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.32]); ylabel(''); % xlabel('Track position (cm)')
% peak5DistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
% set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
% xlabel('Track Position (cm)')

% % Stack all together
% ecSortMap = [tmpMap2; tmpMap3; tmpMap5];
% nUnits = [size(tmpMap5,1) size(tmpMap3,1) size(tmpMap2,1)];
% wflFigLyrs = figure; hold on
% set(gcf,'units','normalized','position',[0.4 0.05 0.25 0.8])
% imagesc(ecSortMap,[prctile(ecSortMap,1,'all'), prctile(ecSortMap,98,'all')]);
% plot([0.5 nBins],[nUnits(3) nUnits(3)],'k--','LineWidth',2)
% plot([0.5 nBins],[sum(nUnits(2:3)) sum(nUnits(2:3))],'k--','LineWidth',2)
% plot([datStruc.rwdBin datStruc.rwdBin],[0 sum(nUnits)],'w--','LineWidth',2)
% colormap("parula")
% cbar = colorbar; clim([0 0.98]);
% xticks(1:10:nBins+1)
% xticklabels(100*datStruc.binedges(1:10:end))
% xlim([0.15 nBins+0.85]); ylim([0.15 sum(nUnits)+0.85])
% ylabel('Unit #'); ylabel(cbar,'FR (Normalized)','FontSize',12,'Rotation',90)
% xlabel('Track Position (cm)')
% set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

% Plot distros together
peakLyrDistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
plot(binpos*100,ec2FieldDistro./sum(ec2FieldDistro),'Color','m','LineWidth',2);
plot(binpos*100,ec3FieldDistro./sum(ec3FieldDistro),'Color','c','LineWidth',2);
plot(binpos*100,ec5FieldDistro./sum(ec5FieldDistro),'Color','k','LineWidth',2);
ylim([0 inf])
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.14])
xlabel('Track Position (cm)')

if saveFlag
    fsave(wflFig2,[sbase '_ec_pilot_waterfall_EC2'],1,1)
    fsave(wflFig3,[sbase '_ec_pilot_waterfall_EC3'],1,1)
    fsave(wflFig5,[sbase '_ec_pilot_waterfall_EC5'],1,1)
    % fsave(wflFigLyrs,[sbase '_ec_pilot_waterfall_EC2-5'],1,1)
    fsave(peakLyrDistroFig,[sbase '_ec_pilot_waterfall_EC2-5_histo'],1,1)
end

%% Putative 1D grid cells
% From Wen et al., 2024
% Z-score spatial FR for each unit
% Compute spectrogram using 16lap / 30m window
% Look for 3 peak structure in PSD
dbnsz = 0.01;
dend = max(sess.pos(sess.lapInclude));
[tmpedges, tmpbinfr] = get_frXpos(root,1,sess,dbnsz,dend);
sfs = size(tmpbinfr,2);
wn = (sfs) * 10;     % 16 laps worth of segments
ovlp = round(wn*0.875);

clear pxxA pxxF pxxN
for i = 1:length(root.good)
    cc = root.good(i);
    % cc = 15;

    [~,binFRA] = get_frXpos(root,cc,sess,dbnsz);
    [~,binFRF] = get_frXpos(rootFrst,cc,sessFrst,dbnsz);
    [~,binFRN] = get_frXpos(rootLast,cc,sessLast,dbnsz);
    
    spFRA = reshape(binFRA',[numel(binFRA),1]);
    spFRA(isnan(spFRA)) = 0;
    spFRA = normalize(spFRA);
    spFRF = reshape(binFRF',[numel(binFRF),1]);
    spFRF(isnan(spFRF)) = 0;
    spFRF = normalize(spFRF);
    spFRN = reshape(binFRN',[numel(binFRN),1]);
    spFRN(isnan(spFRN)) = 0;
    spFRN = normalize(spFRN);

    [pxxA(i,:),f] = xcorr(spFRA,sfs*2); % Compute xcorr over 2 cycles
    [pxxF(i,:)] = xcorr(spFRF,sfs*2);
    [pxxN(i,:)] = xcorr(spFRN,sfs*2);
    
    % Try with spectral power estimate instead 
    % [pxxA,f] = pwelch(spFRA,wn,ovlp,10000,sfs); % f is in m^-1 so place cell will have peak frequency every 1 cycle
    % [pxxF] = pwelch(spFRF,wn,ovlp,10000,sfs);
    % [pxxN] = pwelch(spFRN,wn,ovlp,10000,sfs);
    % 
    % figure; hold on;
    % plot(f,pxx,'k')
    % plot(f,pxxfrst,'g')
    % plot(f,pxxlast,'r')
    % ylim([0 2]); xlim([0 5])
    % ylabel('PSD'); xlabel('spatial freq')
    % title(['Unit ' num2str(cc)])
end

%% Visualize spatial acgs
npxxA = smoothdata(pxxA,2,'movmean',10);
npxxA = normalize(npxxA,2,'range',[0 1]);
npxxF = smoothdata(pxxF,2,'movmean',10);
npxxF = normalize(npxxF,2,'range',[0 1]);
npxxN = smoothdata(pxxN,2,'movmean',10);
npxxN = normalize(npxxN,2,'range',[0 1]);

for i = 1:length(root.good)
    [pks,locs] = findpeaks(npxxA(i,:),'MinPeakProminence',0.1);
    npksA(i) = numel(pks);
    pkDA(i) = mean(diff(locs));
    [pks,locs] = findpeaks(npxxF(i,:),'MinPeakProminence',0.1);
    npksF(i) = numel(pks);
    pkDF(i) = mean(diff(locs));
    [pks,locs] = findpeaks(npxxN(i,:),'MinPeakProminence',0.1);
    npksN(i) = numel(pks);
    pkDN(i) = mean(diff(locs));
end
putGC = npksA >= 5; % # peaks for putative Grid Cells over 2 cycles
[~,inds] = sort(pkDA(putGC));

tmpbns = sfs:1:3*sfs;   % Only plot 1 cycle
allGCF = plot_unitWaterfall(npxxA(putGC,tmpbns),(tmpbns-sfs*2)/100,inds,1,0);
colormap jet; set(gcf,'Position',[.4 .35 .2 .3]); clim([0 0.9])

% fsave(allGCF,[sbase '_GC_acg'])

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

% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%

%% Combine Sessions/Animals
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%

parentDir = "D:\Data\Kelton\analyses\group_analyses"; 
datT = import_xldat(parentDir,"dat_include.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\mEC_Op_Layers';
cd(groupSDir) 

% mInclude = {'KW043','KW052','KW054','KW055','KW093','KW098'}; % For behavior
mInclude = {'KW043','KW093','KW098'}; %For ephys
dInclude = 1:3;
% useInds = datT.include == 1;
for i = 1:height(datT)
    useInds(i) = logical(sum(strcmp(datT.mouse(i),mInclude))) & logical(sum(ismember(dInclude,datT.session(i)))); 
end
datT(~useInds,:) = [];  %Clean excluded sessions

saveFlag = 1;
sbase = 'mEC_pilot_';
fname = [sbase 'physDat1'];

dbnsz = 0.05;
histoBnsz = 5;
binpos = 0.025:dbnsz:1.825;
binedges = 0:0.05:1.85;
wlen = 150;
r1pos = 0.9;    % 10 cm
vColors2 = [0.5 0.5 1; 0.75 0.75 1];
r1posInd = find(binpos > r1pos,1);

%% Combine EC Behavior data
combine_bhvrDat(datT,fname,groupSDir);

%% Combine EC Phys data
combine_ephysDat(datT,fname,groupSDir);

%% Load old data
% load('D:\Data\Kelton\analyses\group_analyses\dp2_progreport\dp2_progReport_data.mat');
load(fullfile(groupSDir,fname))

nSess = size(recID,1);
nMice = numel(unique(recID(:,1)));

%% Analyze group behavior data

for i = 1:nSess
    uSpVel(i,:) = mean(bvDat(i).spVelMap,1,'omitnan');
    uSpLck(i,:) = mean(bvDat(i).lckRMap,1,'omitnan');
    uLckDI(i,:) = bvDat(i).ulckDI;
end
for i = 1:length(dInclude)
    dayStruc(i).uSpVel = mean(uSpVel(recID(:,2) == i,:));
    dayStruc(i).uSpLck = mean(uSpLck(recID(:,2) == i,:));
    [dayStruc(i).ciup,dayStruc(i).cidn] = get_CI(uSpLck(recID(:,2) == i,:));
    dayStruc(i).uLckDI = mean(uLckDI(recID(:,2) == i,:));
end

% % RM Anova (Apparently 2024b enbables functionality for multcompare on ranova 
% ldiXday = reshape(uLckDI,[length(dInclude) nMice])';
% ldiT = table(unique(recID(:,1)),ldiXday(:,1),ldiXday(:,2),ldiXday(:,3),'VariableNames',{'mouse','d1','d2','d3'});
% rm = fitrm(ldiT,'d1-d3 ~ 1','WithinDesign',[1:3]);
% rmanova = ranova(rm,'WithinModel',[1:3]');
% multcompare(rma,[1:3]')

% Plot spatially binned lick rate
vColors = [0 0 0; 0 0 1; 0.75 0.75 1];
binedges = 0.015:0.03:1.835;
lckPlotXDay = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.5 0.22 0.23])
for i = 1:length(dInclude)
    plot_CIs(binedges,dayStruc(i).ciup,dayStruc(i).cidn,vColors(i,:))
    plot(binedges,dayStruc(i).uSpLck,'Color',vColors(i,:))
end
plot([r1pos r1pos],[0 10],'--','Color',vColors(2,:))
xlabel('Position (cm)'); xlim([0 1.85])
ylabel('Lick Rate (Hz)')
set(gca,'FontSize',16,'FontName','Arial')

% Plot LDI by day
ldiPlotXDay = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.1 0.27])
b = bar(dInclude,vertcat(dayStruc.uLckDI),'FaceColor','flat','HandleVisibility','off');
b.CData = vColors;
for i = 1:length(dInclude)
    errorbar(i,dayStruc(i).uLckDI,std(vertcat(uLckDI(recID(:,2) == i,:)))/sqrt(nMice),'k.','HandleVisibility','off')
end
plot(reshape(uLckDI,[length(dInclude) nMice]),'k-o')
xlim([0.5 3.5])
xticks(1:3)
ylabel('Lick Selectivity')
% legend({'Mouse'},'location','southeast')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag 
    fsave(lckPlotXDay,[sbase 'bvr_lckXday'],1,1)
    fsave(ldiPlotXDay,[sbase 'bvr_ldiBar'],1,1)
end
%% 
% recID = [];     % [mouseID, recDay, recUnitID, recLayer]
% useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
% lcDat = [];     % [frstHalf.si_p, frstHalf.si, frstHalf.pkFR, frstHalf.pkLoc, lastHalf.si_p, lastHalf.si, lastHalf.pkFR, lastHalf.pkLoc]
% lcMap = [];     % [frstHalf.posfr, lastHalf.posfr];
% siStr = [];     % Struc containing binned SI data per 10 laps
% % vlDat = [];     % [frstHalf sig., frstHalf slope, frstHalf R, lastHalf sig., lastHalf slope, lastHalf R]
% % thDat = [];     % [frstHalf.p, frstHalf.mrl, frstHalf.ang, lastHalf.p, lastHalf.mrl, lastHalf.ang]
% % thMap = [];     % [frstHalf.thetafr, lastHalf.thetafr];
% % bvDat = [];     % [frstHalf.lckDI, frstHalf.preRZV, lastHalf.lckDI, lastHalf.preRZV
% 
% clear ps stats
% 
% ct    = 1;
% 
% for i = 1:height(datT)
%     if datT.sess_type(i) ~= 6
%         continue
%     end
% 
%     % === Load data ===
%     cd(datT.fpath{i})
% 
%     rootfile = dir("*_root.mat");
%     load(rootfile.name)
%     sessfile = dir("*_session.mat");
%     load(sessfile.name)
%     epochfile = dir("*_ec_pilot_dat.mat");
%     load(epochfile.name)
%     disp(root.name)
% 
%     useUnits = root.info.fr(root.goodind) > 0.1 & root.info.uType(root.goodind);
%     nCCs = length(root.good);
% 
%     % === Concatenate recording data ===
%     useCC = logical([useCC; useUnits]);
%     recID = [recID; str2num(datT.mouse{i}(end-2:end))*ones(nCCs,1), datT.session(i)*ones(nCCs,1), root.good, datStruc.lyrID];
% 
%     % === Concatenate SI and Peak data ===
%     lcDat = [lcDat; datStruc.sigSI, datStruc.trueSI, datStruc.truePk, datStruc.trueLc];
%     lcMap = [lcMap; datStruc.posfr];
% 
%     ct = ct + 1;
% end
% 
% cd(groupSDir)
% 
% save([sbase 'data'],'recID','useCC','lcDat','lcMap','siStr')

%% Waterfall across days using all mice

for i = 1:max(celID(:,2))
    useEC2 = useCC & celID(:,4) == 2 & lcDat(:,1) < 0.05 & celID(:,2) == i;
    useEC3 = useCC & celID(:,4) == 3 & lcDat(:,1) < 0.05 & celID(:,2) == i;
    useEC5 = useCC & celID(:,4) == 5 & lcDat(:,1) < 0.05 & celID(:,2) == i;

    [wflFig2,tmpMap2,sortEC2] = plot_unitWaterfall(lcMap(useEC2,:),binedges,0,1,0);
    set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.35])
    ec2FieldDistro = histcounts(binpos(lcDat(useEC2,4)),binedges);
    [wflFig3,tmpMap3,sortEC3] = plot_unitWaterfall(lcMap(useEC3,:),binedges,0,1,0);
    set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.35])
    ec3FieldDistro = histcounts(binpos(lcDat(useEC3,4)),binedges);
    [wflFig5,tmpMap5,sortEC5] = plot_unitWaterfall(lcMap(useEC5,:),binedges,0,1,0);
    set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.35])
    ec5FieldDistro = histcounts(binpos(lcDat(useEC5,4)),binedges);
    
    pkDistroXLyr = plotDistroHisto([ec2FieldDistro; ec3FieldDistro; ec5FieldDistro],binpos,r1pos);
    set(gcf,'units','normalized','position',[0.4 0.35 0.21 0.14])
    xlabel('Track Position (cm)')

    if saveFlag
        fsave(wflFig2,[sbase 'wfl_D' num2str(i) '_EC2'],1,1)
        fsave(wflFig3,[sbase 'wfl_D' num2str(i) '_EC3'],1,1)
        fsave(wflFig5,[sbase 'wfl_D' num2str(i) '_EC5'],1,1)
        fsave(pkDistroXLyr,[sbase 'wfl_Distro_D' num2str(i)],1,1)
        close all
    end
end

%% Plot Ratio of sig. spatial units by day, by layer
% Organized as mouse, layer, day
mID = unique(celID(:,1));
for i = 1:max(celID(:,2))
    for j = 1:length(mID)
        nSI(j,1,i) = sum(useCC & celID(:,4) == 2 & lcDat(:,1) < 0.05 & celID(:,2) == i & celID(:,1) == mID(j));
        nSI(j,2,i) = sum(useCC & celID(:,4) == 3 & lcDat(:,1) < 0.05 & celID(:,2) == i & celID(:,1) == mID(j));
        nSI(j,3,i) = sum(useCC & celID(:,4) == 5 & lcDat(:,1) < 0.05 & celID(:,2) == i & celID(:,1) == mID(j));
        nCC(j,1,i) = sum(useCC & celID(:,4) == 2 & celID(:,2) == i & celID(:,1) == mID(j));
        nCC(j,2,i) = sum(useCC & celID(:,4) == 3 & celID(:,2) == i & celID(:,1) == mID(j));
        nCC(j,3,i) = sum(useCC & celID(:,4) == 5 & celID(:,2) == i & celID(:,1) == mID(j));
    end
end

siRatio = nSI ./ nCC;

% Plot siRatio by day

siRatioEC2 = plot_datXday(squeeze(siRatio(:,1,:)),[1,0,1]);
siRatioEC3 = plot_datXday(squeeze(siRatio(:,2,:)),[0,1,1]);
siRatioEC5 = plot_datXday(squeeze(siRatio(:,3,:)),[.9 .9 .9]);

if saveFlag
    fsave(siRatioEC2,[sbase, 'siRatio_EC2'],1,1)
    fsave(siRatioEC3,[sbase, 'siRatio_EC3'],1,1)
    fsave(siRatioEC5,[sbase, 'siRatio_EC5'],1,1)
end

%% Outdated waterfall plotting

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
nBins = length(binpos);

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

function [fhandle] = plotDistroHisto(distro,binpos,rzPos)
cols = {'m','c','k'};

fhandle = figure; hold on
for i = 1:size(distro,1)
plot(binpos*100,distro(i,:)./sum(distro(i,:)),'Color',cols{i},'LineWidth',2);
end
plot([rzPos rzPos]*100,[0 0.15],'k--','HandleVisibility','off')
% ylim([0 max(distro./sum(distro),[],'all')])
xlim([binpos(1)*100-1 binpos(end)*100+1])
ylabel('P(field peak)')
set(gca,'FontSize',12,'FontName','Arial')
end

function [fhandle] = plot_datXday(dat,vCols)
nDays = size(dat,2);
nMice = size(dat,1);

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.35])
b = bar(1:nDays,mean(dat),'FaceColor','flat','HandleVisibility','off');
b.CData = [vCols/4;vCols/2;vCols];
errorbar(1:nDays,mean(dat),std(dat)/sqrt(nMice),'k.','HandleVisibility','off')
plot(dat','k-o')
xlim([0.5 nDays+0.5])
xticks(1:nDays);
ylabel('P(Spatial)')
set(gca,'FontSize',12,'FontName','Arial')
end
