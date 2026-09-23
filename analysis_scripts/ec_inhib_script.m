spath = 'D:\Data\Kelton\analyses\KW109\KW109_07222026_rec_D2_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW111\KW111_08072026_rec_D2_RLat2';
% spath = 'D:\Data\Kelton\analyses\KW112\KW112_08182026_rec_D2_RMed1';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_dat.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

rwdShift = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
rwdShift = sess.rwdTrials(rwdShift);
nUnits = length(root.good);
saveFlag = 1;

dbnsz = 0.05;
binpos = dbnsz/2:dbnsz:1.85-dbnsz/2;

%% Calculate P(Rwd) at Nov location in Fam

pRwdFam_atNov = zeros(length(sessFrst.valTrials),1);
rz_nov = [1.3 1.6];

for i = 1:length(sessFrst.valTrials)
    tmpstt = sessFrst.lapstt(sessFrst.valTrials(i));
    tmpend = sessFrst.lapend(sessFrst.valTrials(i));
    tmplck = sessFrst.lckind(sessFrst.lckind > tmpstt & sessFrst.lckind < tmpend);
    tmpcts = histcounts(sessFrst.pos(tmplck),rz_nov);
    pRwdFam_atNov(i) = tmpcts > 0;
end

disp(['P(Rwd) in Fam at Nov location ' num2str(mean(pRwdFam_atNov))]);
disp(['P(Rwd) in Nov at Nov location ' num2str(length(sessLast.rwdTrials) / length(sessLast.valTrials))]);

%% Plot opto graphs individually

cd(spath)
mkdir('optoPlots')
cd('optoPlots')

sess.optoind = sess.optoind(1:205);
for i = 1:length(root.good)
    cc = root.good(i);

    [~,~,tmpOptoF] = plot_frXopto(root,cc,sess);

    saveas(tmpOptoF,['unit',num2str(cc),'_opto'],'png')
    close(tmpOptoF)
end

cd(spath)

%% Recalculate opto after plotting it separately

[~,sess.optoind]= findpeaks(double(sess.opto > 2));
sess.opto = sess.opto > 2;
sess.optoUpInd = sess.ind(sess.opto);

%% Find opto laps
sess.optolapinds = zeros(size(sess.ts))';
sess.optolap = logical(zeros(size(sess.valTrials)));

for i = 1:length(sess.valTrials)
    tmpopto = sess.optoUpInd(sess.optoUpInd > sess.lapstt(i) & sess.optoUpInd < sess.lapend(i));
    optobins = histcounts(sess.pos(tmpopto),0:0.05:1.85) > 1;
    if sum(optobins) > length(optobins) / 2
        sess.optolapinds(sess.lapstt(i):sess.lapend(i)) = true;
        sess.optolap(i) = true;
    else
        sess.optolap(i) = false;
    end
end

tOpto = sum(sess.optolapinds) ./ sess.samprate;
tBase = sess.ts(end) - tOpto;

%% Calculate burst, index rate, and length in opto laps vs non opto laps
bstNum = zeros(length(root.good),2);
bstLen = zeros(length(root.good),2);
spkNum = zeros(length(root.good),2);

for i = 1:length(root.good)
    cc = root.good(i);
    bstIndex(i) = get_burstIndex(root,sess,cc);

    tmpSpk = root.tsb(root.cl == cc);
    useBst = find(root.burst_cl == cc);
    tmpBst = root.burst_tsb(useBst);
    baseBsts = find(sess.optolapinds(tmpBst) == 0);
    optoBsts = find(sess.optolapinds(tmpBst) == 1);
    % baseBsts = tmpbst(sess.optolapinds(tmpbst) == 0);
    % optoBsts = tmpbst(sess.optolapinds(tmpbst) == 1);

    bstNum(i,:) = [numel(baseBsts) numel(optoBsts)];
    bstLen(i,:) = [mean(root.burst_len(useBst(baseBsts))) mean(root.burst_len(useBst(optoBsts)))];
    spkNum(i,:) = [sum(sess.optolapinds(tmpSpk)) sum(~sess.optolapinds(tmpSpk))];
end

bstRate = [bstNum(:,1)./tBase bstNum(:,1)./tOpto];
spkRate = [spkNum(:,1)./tBase spkNum(:,1)./tOpto];

%% Compare Sub pyr burst rate
useUnits = root.info.lyrID(root.goodind) == 1 & root.info.uType(root.goodind) & bstIndex' > 0.1;
% useINs = root.info.lyrID(root.goodind) == 1 & ~root.info.uType(root.goodind);

[~,ps.bstRate] = ttest(bstRate(useUnits,1),bstRate(useUnits,2));
[~,ps.bstLen] = ttest(bstLen(useUnits,1),bstLen(useUnits,2));
[~,ps.spkRate] = ttest(spkRate(useUnits,1),spkRate(useUnits,2));
[~,ps.spkNormBstRate] = ttest(bstRate(useUnits,1)./spkRate(useUnits,1),bstRate(useUnits,2)./spkRate(useUnits,2));

% % Visualize IN vs Pyr bursts
% figure; hold on;
% plot(ones(sum(useINs),1),bstIndex(useINs),'r.')
% plot(2*ones(sum(useUnits),1),bstIndex(useUnits),'b.')

bstRateF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
bstRateF = fixRatio(bstRateF);
bar(mean(bstRate(useUnits,:),"omitmissing"));
plot(bstRate(useUnits,:)','k-o')
xticks([1 2]); xlim([.5 2.5]);
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(bstRateF,'Burst Rate',ps.bstRate)

bstLenF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
bstLenF = fixRatio(bstLenF);
bar(mean(bstLen(useUnits,:),"omitmissing"));
plot(bstLen(useUnits,:)','k-o')
xticks([1 2]); xlim([.5 2.5]); ylim([3 5])
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(bstLenF,'Burst Length',ps.bstLen)

spkRateF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
spkRateF = fixRatio(spkRateF);
bar(mean(spkRate(useUnits,:),"omitmissing"));
plot(spkRate(useUnits,:)','k-o')
xticks([1 2]); xlim([.5 2.5]);
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(spkRateF,'Spike Rate',ps.spkRate)
%%
spkNormBstRateF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
spkNormBstRateF = fixRatio(spkNormBstRateF);
bar(mean(bstRate(useUnits,:)./spkRate(useUnits,:),"omitmissing"));
plot((bstRate(useUnits,:)./spkRate(useUnits,:))','k-o')
xticks([1 2]); xlim([.5 2.5]);
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(spkNormBstRateF,'Burst Rate / Spike Rate',ps.spkNormBstRate)
%%
if saveFlag 
    fsave(bstRateF,[root.name '_burstRateComp'],1,0)
    fsave(bstLenF,[root.name '_burstLengthComp'],1,0)
    fsave(spkRateF,[root.name '_spikeRateComp'],1,0)
end

%% Reward Shift waterfall analyses - whole session

lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = lyrUnits & hiFRUnits & root.info.uType(root.goodind);
frstHalf.useUnits = useUnits; lastHalf.useUnits = useUnits;
siUnits = useUnits & (frstHalf.sigSI <= 0.05 | lastHalf.sigSI <= 0.05);
frstSIUnits = useUnits & (frstHalf.sigSI <= 0.05);
lastSIUnits = useUnits & (lastHalf.sigSI <= 0.05);
bothSIUnits = useUnits & (lastHalf.sigSI <= 0.05 & frstHalf.sigSI <= 0.05);

preRZBin = find(frstHalf.binedges > 0.4,1);
pstRZBin = find(frstHalf.binedges > 1.3,1);

famWflF = plot_unitWaterfall(frstHalf.posfr(frstSIUnits,:),frstHalf.binedges,0,1,0);
plot([preRZBin preRZBin],[0 sum(frstSIUnits)],'k--')

novWflF = plot_unitWaterfall(lastHalf.posfr(lastSIUnits,:),frstHalf.binedges,0,1,0);
plot([pstRZBin pstRZBin],[0 sum(lastSIUnits)],'r--')

% Both SI units
[famWflBothF,~, famSort] = plot_unitWaterfall(frstHalf.posfr(bothSIUnits,:),frstHalf.binedges,0,1,0);
plot([preRZBin preRZBin],[0 sum(bothSIUnits)],'k--')

novWflBothF = plot_unitWaterfall(lastHalf.posfr(bothSIUnits,:),frstHalf.binedges,famSort,1,0);
plot([pstRZBin pstRZBin],[0 sum(bothSIUnits)],'r--')

if saveFlag
    fsave(famWflF,[root.name '_lc_wfl_frstSI_fam'],1,0);
    fsave(novWflF,[root.name '_lc_wfl_lastSI_nov'],1,0);
    fsave(famWflBothF,[root.name '_lc_wfl_bothSI_fam'],1,0);
    fsave(novWflBothF,[root.name '_lc_wfl_bothSI_nov'],1,0);
end

%% Reward Shift waterfall analyses - opto laps vs non session

frstHalf.optoLaps = sess.optolap(1:rwdShift-1); % Accounts for 1st trial non valid
lastHalf.optoLaps = sess.optolap(rwdShift:end);
uPosFrBaseFam = squeeze(mean(frstHalf.frMap(~frstHalf.optoLaps,:,:)))';
uPosFrOptoFam = squeeze(mean(frstHalf.frMap(frstHalf.optoLaps,:,:)))';
uPosFrBaseNov = squeeze(mean(lastHalf.frMap(~lastHalf.optoLaps,:,:)))';
uPosFrOptoNov = squeeze(mean(lastHalf.frMap(lastHalf.optoLaps,:,:)))';

famWflF_base = plot_unitWaterfall(uPosFrBaseFam(bothSIUnits,:),frstHalf.binedges,0,1,0);
plot([preRZBin preRZBin],[0 sum(frstSIUnits)],'k--')
title('Fam. Baseline')
famWflF_opto = plot_unitWaterfall(uPosFrOptoFam(bothSIUnits,:),frstHalf.binedges,0,1,0);
plot([preRZBin preRZBin],[0 sum(frstSIUnits)],'k--')
title('Fam. Opto')

novWflF_base = plot_unitWaterfall(uPosFrBaseNov(bothSIUnits,:),frstHalf.binedges,0,1,0);
plot([pstRZBin pstRZBin],[0 sum(frstSIUnits)],'r--')
title('Nov. Baseline')
novWflF_opto = plot_unitWaterfall(uPosFrOptoNov(bothSIUnits,:),frstHalf.binedges,0,1,0);
plot([pstRZBin pstRZBin],[0 sum(frstSIUnits)],'r--')
title('Nov. Opto')

if saveFlag
    fsave(famWflF_opto,[root.name '_lc_wfl_bothSI_F_Opto'],1,0);
    fsave(famWflF_base,[root.name '_lc_wfl_bothSI_F_Base'],1,0);
    fsave(novWflF_base,[root.name '_lc_wfl_bothSI_N_Opto'],1,0);
    fsave(novWflF_opto,[root.name '_lc_wfl_bothSI_N_Base'],1,0);
end

%% PVCs 
posNormFamOff = (uPosFrBaseFam(bothSIUnits,:)' ./ max(uPosFrBaseFam(bothSIUnits,:)'))';
posNormFamStm = (uPosFrOptoFam(bothSIUnits,:)' ./ max(uPosFrOptoFam(bothSIUnits,:)'))';
posNormNovOff = (uPosFrBaseNov(bothSIUnits,:)' ./ max(uPosFrBaseNov(bothSIUnits,:)'))';
posNormNovStm = (uPosFrOptoNov(bothSIUnits,:)' ./ max(uPosFrOptoNov(bothSIUnits,:)'))';

pvFamOffStm = corr(posNormFamOff,posNormFamStm);   % units x sp. bins
pvFamOffStmF = plot_pvcorr(pvFamOffStm);
plot([preRZBin preRZBin],[0 sum(bothSIUnits)],'w--')
title('Fam. Baseline vs Opto')
pvNovOffStm = corr(posNormNovOff,posNormNovStm);   % units x sp. bins
pvNovOffStmF = plot_pvcorr(pvNovOffStm);
plot([pstRZBin pstRZBin],[0 sum(bothSIUnits)],'r--')
title('Nov. Baseline vs Opto')

pvFamNovOff = corr(posNormFamOff,posNormNovOff);   % units x sp. bins
pvFamNovOffF = plot_pvcorr(pvFamNovOff);
plot([preRZBin preRZBin],[0 sum(bothSIUnits)],'w--')
title('Baseline Fam v Nov')
pvFamNovStm = corr(posNormFamStm,posNormNovStm);   % units x sp. bins
pvFamNovStmF = plot_pvcorr(pvFamNovStm);
plot([preRZBin preRZBin],[0 sum(bothSIUnits)],'w--')
title('Opto Fam v Nov')

if saveFlag
    fsave(pvFamOffStmF,[root.name '_lc_pvc_bothSI_F_BaseVOpto'],1,0);
    fsave(pvNovOffStmF,[root.name '_lc_pvc_bothSI_N_BaseVOpto'],1,0);
    fsave(pvFamNovOffF,[root.name '_lc_pvc_bothSI_Base_FamVNov'],1,0);
    fsave(pvFamNovStmF,[root.name '_lc_pvc_bothSI_Opto_FamVNov'],1,0);
end

%%
% Make off-diagonal matrix
idMat = logical(eye(nBins));

% PVC fam off vs on for visualization
posNormPreOff = (squeeze(mean(frstHalf.frMap(~sess.optolap(1:rwdShift-1),:,bothSIUnits),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(~sess.optolap(1:rwdShift-1),:,bothSIUnits),1,"omitnan")),[],1))';
posNormPreStm = (squeeze(mean(frstHalf.frMap(sess.optolap(1:rwdShift-1),:,bothSIUnits),1,"omitnan")) ./ max(squeeze(mean(frstHalf.frMap(sess.optolap(1:rwdShift-1),:,bothSIUnits),1,"omitnan")),[],1))';
pvPreOffStm = corr(posNormPreOff,posNormPreStm);   % units x sp. bins
pvPreOffStmCompF = plot_pvcorr(pvPreOffStm);
xlabel('Position (cm) Familiar Off'); ylabel('Position (cm) Familiar Stim');

% PVC nov off vs on for visualization
posNormPstOff = (squeeze(mean(lastHalf.frMap(~sess.optolap(rwdShift:end),:,bothSIUnits),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(~sess.optolap(rwdShift:end),:,bothSIUnits),1,"omitnan")),[],1))';
posNormPstStm = (squeeze(mean(lastHalf.frMap(sess.optolap(rwdShift:end),:,bothSIUnits),1,"omitnan")) ./ max(squeeze(mean(lastHalf.frMap(sess.optolap(rwdShift:end),:,bothSIUnits),1,"omitnan")),[],1))';
pvPstOffStm = corr(posNormPstOff,posNormPstStm);   % units x sp. bins
pvPstOffStmCompF = plot_pvcorr(pvPstOffStm);
xlabel('Position (cm) Novel Off'); ylabel('Position (cm) Novel Stim');

if saveFlag
    fsave(pvPreOffStmCompF,[root.name, '_pvCorrPreOffStm'],1,0);
    fsave(pvPstOffStmCompF,[root.name, '_pvCorrPstOffStm'],1,0);
end
