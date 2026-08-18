spath = 'D:\Data\Kelton\analyses\KW109\KW109_07222026_rec_D2_RLat1';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_dat.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

try
    rwdShift = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
    if isfield(sess,'valTrials')
        rwdShift = sess.valTrials(rwdShift);
    end
catch
end

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
bar(mean(bstRate(useUnits,:),"omitmissing"));
plot(bstRate(useUnits,:)','k-o')
xticks([1 2]); xlim([.5 2.5]);
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(bstRateF,'Burst Rate',ps.bstRate)

bstLenF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
bar(mean(bstLen(useUnits,:),"omitmissing"));
plot(bstLen(useUnits,:)','k-o')
xticks([1 2]); xlim([.5 2.5]); ylim([3 5])
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(bstLenF,'Burst Length',ps.bstLen)

spkRateF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
bar(mean(spkRate(useUnits,:),"omitmissing"));
plot(spkRate(useUnits,:)','k-o')
xticks([1 2]); xlim([.5 2.5]);
xticklabels({'Baseline','Opto'})
set(gca,'FontSize',16,'FontName','Arial')
text2bar(spkRateF,'Spike Rate',ps.spkRate)
%%
spkNormBstRateF = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.14 0.4])
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

%% PV On stim vs off, pre and post
frstHalf.sigSI = sum(frstHalf.shufSI >= frstHalf.trueSI,2) / nShufs;
lastHalf.sigSI = sum(lastHalf.shufSI >= lastHalf.trueSI,2) / nShufs;

nBins = length(binpos);
lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = lyrUnits & hiFRUnits & root.info.uType(root.goodind);
bothSIUnits = useUnits & (lastHalf.sigSI <= 0.05 & frstHalf.sigSI <= 0.05);
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
