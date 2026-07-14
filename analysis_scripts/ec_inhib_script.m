spath = 'D:\Data\Kelton\analyses\KW106\KW106_06252026_rec_D2_RLat1';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_dat.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

try
    rwdShift = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
    if isfield(sess,'rwdTrials')
        rwdShift = sess.rwdTrials(rwdShift);
    end
catch
end

nUnits = length(root.good);
saveFlag = 1;

%% Recalculate opto in case of licks contamination

[~,sess.optoind]= findpeaks(double(sess.opto > 2));
sess.opto = sess.opto > 2;

%% Find opto laps
sess.optolapinds = zeros(size(sess.ts))';

for i = 1:length(sess.valTrials)
    tmpopto = sess.optoind(sess.optoind > sess.lapstt(i) & sess.optoind < sess.lapend(i));
    optobins = histcounts(sess.pos(tmpopto),0:0.05:1.85) > 1;
    if sum(optobins) > length(optobins) / 2
        sess.optolapinds(sess.lapstt(i):sess.lapend(i)) = true;
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

if saveFlag
    fsave(bstRateF,[root.name '_burstRateComp'])
    fsave(bstLenF,[root.name '_burstLengthComp'])
    fsave(spkRateF,[root.name '_spikeRateComp'])
end
