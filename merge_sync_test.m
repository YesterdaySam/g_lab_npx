%% Merge data streams

sess = load("KW001_05302024_session.mat");

root = load("KW001_05302024_root.mat");

sync = load("KW001_05302024_sync.mat");
sync.ts = [1:length(digArray)]/2500;

%%

% recstart = find(sess.slx < 1,1);
recend = find(sess.slx < 1,1,'last');
recstart = sess.ts(recend) - str2double(meta.fileTimeSecs);
recstart = find(sess.ts > recstart,1);

sess2 = epochsess(sess,[recstart recend]);

%% unaligned
figure; plot([1:length(sync.digArray)]/2500,sync.digArray)
hold on; 
plot([1:length(sess2.slx)]/2000,0.1+[sess2.slx>1])

%% aligned
syncend = sync.ts(find(sync.digArray < 0.5, 1, 'last'));
[~,locs] = findpeaks(double(sess2.slx > 1));
reclastpulse = sess.ts(locs(end)-1);

offset = diff([syncend reclastpulse]);

figure; plot([1:length(sync.digArray)]/2500 + offset,sync.digArray)
hold on; 
plot([1:length(sess2.slx)]/2000,0.1+[sess2.slx>1])

%% Plot velocity-binned firing rate

cc = 50;

ccspk = root.ts(root.cl == cc) + offset;
ccspk = ccspk(ccspk > 0);

vbnsz = 0.02; %m/s

vedges = 0:vbnsz:max(sess2.vel);

for i = 1:length(ccspk)
    spkvel(i) = sess2.vel(find(sess2.ts > ccspk(i),1));
end

vspk = histcounts(spkvel,vedges);
vocc = histcounts(sess2.vel,vedges)/sess2.samprate;

vfr = vspk./vocc;
mdl = fitlm(vedges(1:end-1), vfr);
ys = predict(mdl,[min(vedges); max(vedges)]);

figure; hold on
plot([vedges(1:end-1)]*100,vfr, 'ko','MarkerFaceColor','k')
plot([min(vedges); max(vedges)]*100,ys,'r','LineWidth',2)
xlabel('Velocity (cm/s)'); ylabel('Firing Rate')
title(['Unit ' num2str(cc)])

saveas(gcf,['velXFR_unit' num2str(cc)],'png')

%% Merge data streams #2 - flipper pulse

load("KW004_06272024_session.mat");
sess.slx = double(sess.slx > 0.5);

load("KW004_06272024_root.mat");

sync = load("KW004_06272024_sync.mat");
sync.digArray = double(sync.digArray);
sync.ts = [1:length(digArray)]/2500;

if sess.samprate ~= 2500
    sess.dsslx = downsample(sess.slx,sess.samprate/2500);
    sess.dsts  = downsample(sess.ts,sess.samprate/2500);
    [sync_ccg, sync_lags] = xcorr(sync.digArray, sess.dsslx);
else
    [sync_ccg, sync_lags] = xcorr(sync.digArray, sess.slx);
end

pk_lag = sync_lags(sync_ccg == max(sync_ccg));  %Index to offset sess.ts

%% TShift sess to match sync
%Find sync edge inds and ts
[~, edges_sync] = findpeaks(sync.digArray);
edges_sync_ts = [sync.ts(1) sync.ts(edges_sync)];   %Start at sync.ts wall time (not a proper down/up flip)
tmp = diff(edges_sync_ts);

%Find downsampled sess edge inds and ts
[~, edges_sess] = findpeaks(sess.dsslx);
% edges_sess = edges_sess+pk_lag;
align_start = edges_sess+pk_lag > 0;
edges_sess = edges_sess(align_start);
edges_sess_dsts = sess.dsts(edges_sess);

% [~, edges_sess] = findpeaks(sess.slx);
% align_start = edges_sess+pk_lag > 0;
% edges_sess = edges_sess(align_start);
% edges_sess_ts = sess.ts(edges_sess);

%Find full rate sess edge inds and ts
[~, edges_sess] = findpeaks(sess.slx);
edges_sess_ts = sess.ts(edges_sess);
edges_sess_ts = edges_sess_ts(edges_sess_ts - sess.dsts(abs(pk_lag)) > 0);
edges_sess_ts = [edges_sess_ts(1)-tmp, edges_sess_ts];

%% example alignment of all spike times to behavior
% SA = EA + RA*(SB - EB)/RB
% SA is scaled, bounded, timestamp of event from B mapped to stream A
% SB is timestamp of event from B
% E is TS of Nearest Edge in A or B
% R is sampling rate of stream A or B

allts = root.ts;
root.tssync = allts;

for i = 1:length(allts)
    eb = edges_sync_ts(find(edges_sync_ts < allts(i),1,'first'));
    ea = edges_sess_ts(find(edges_sess_ts > eb,1,'first'));
    try
        root.tssync(i) = ea + 2500*(allts(i) - eb)/2500;
    catch
        root.tssync(i) = allts(i);
    end
end

%% Visualize
figure; hold on; plot(sync.ts,sync.digArray)
plot(sess.dsts+pk_lag/2500,sess.dsslx/2)
plot(sess.ts+pk_lag/2500,sess.slx/3)

%% better spkind method from CMBHome

% split spikes into bins with edges of sess.ts
counts = histcounts(tprime,sess.ts);   %after aligning all spike ts use root.ts
inds = [];

% find index of all bins >i spikes and concatenate to existing indices
for i = 0:max(counts)
    inds = cat(2,inds,find(counts>i));
end
% sort inds in case of repeats
inds = sort(inds);
inds = inds(:);

root.tsb = inds;

%% Plot velocity-binned firing rate

cc = 2;

vbnsz = 0.02; %m/s

vedges = 0:vbnsz:max(sess.velshft);

for i = 1:length(tprime)
    spkvel(i) = sess.velshft(find(sess.ts > tprime(i),1));
end

vspk = histcounts(spkvel,vedges);
vocc = histcounts(sess.velshft,vedges)/sess.samprate;

vfr = vspk./vocc;
mdl = fitlm(vedges(1:end-1), vfr);
ys = predict(mdl,[min(vedges); max(vedges)]);

figure; hold on
plot([vedges(1:end-1)]*100,vfr, 'ko','MarkerFaceColor','k')
plot([min(vedges); max(vedges)]*100,ys,'r','LineWidth',2)
xlabel('Velocity (cm/s)'); ylabel('Firing Rate')
title(['Unit ' num2str(cc)])

% saveas(gcf,['velXFR_unit' num2str(cc)],'png')
