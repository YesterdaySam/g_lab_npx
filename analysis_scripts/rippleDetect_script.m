% Dummy LFP / Ripple analysis
% Based off Jadhav lab methods and Buzsaki lab code
% https://www.cell.com/neuron/fulltext/S0896-6273(19)30785-8
% https://github.com/buzsakilab/buzcode/blob/master/analysis/SharpWaveRipples/bz_FindRipples.m
% 
% Filter 150-250 Hz
% Get envelope of signal with Hilbert Transform
% Detect ripple peaks at >3SD of mean on that electrode
% Start/stop as times when envelope exceeds mean near each event
% Optionally discard excessively short/long ripples
% Optionally discard events during runs >4cm/s

saveFlag = 1;

% chan = find(root.lfpinfo.lfpch == 336);
chan = root.uPSDMax(2,3);

root.ripples = get_ripples(root,chan,sess,3,5,[15 250]);

rawlf = root.lfp(chan,:);
riplf = bandpass(rawlf, [150, 250], root.fs_lfp);

%% Optional plotting

% Raw, filtered, and ripple event times over all session
tmplfpfig = figure; hold on    
set(gcf,'units','normalized','position',[0.4 0.35 0.5 0.39])
plot(sess.ts(root.lfp_tsb), rawlf,'k')
plot(sess.ts(root.lfp_tsb), 2+riplf,'r')
plot(sess.ts(root.ripples(:,2)), 2.5 * ones(size(root.ripples,1),1), 'b*')
xlabel('Time (s)')
legend({'Raw','150-250Hz','Ripple Pks'})
set(gca,'FontSize',12,'FontName','Arial')

% Ripple durations
ripStts = sess.ts(root.ripples(:,1));
ripStps = sess.ts(root.ripples(:,3));
ripDurs = ripStps - ripStts;
tmpRipDfig = figure;
histogram(ripDurs*1000)
xlim([0 250])
xlabel('Ripple Duration (ms)')
ylabel('Counts')
set(gca,'FontSize',12,'FontName','Arial')

% Ripple peak power
tmpRipPfig = figure;
histogram(root.ripples(:,4))
xlabel('Ripple Peak Power Amplitude')
xlim([0 0.25])
ylabel('Counts')
set(gca,'FontSize',12,'FontName','Arial')

% Inter-Ripple Interval figures
IRIs = diff(root.ripples(:,1)) ./ root.fs_lfp;

tmpIRIfig = figure;
histogram(IRIs, 0:1:max(IRIs))
% boxplot(IRIs)
xlabel('Inter-Ripple-Interval (s)')
ylabel('Counts')
set(gca,'FontSize',12,'FontName','Arial')

tmpIRIfigZoom = figure;
histogram(IRIs, 0:0.1:max(IRIs))
xlabel('Inter-Ripple-Interval (s)')
ylabel('Counts')
xlim([0 5])
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    sbase = ['chan' num2str(root.lfpinfo.lfpch(chan)) '_sh' num2str(root.lfpinfo.lfpShank(chan))];
    saveas(tmplfpfig,[sbase '_LFP.png'])
    saveas(tmpRipDfig,[sbase '_rippleDur.png'])
    saveas(tmpRipPfig,[sbase '_ripplePwr.png'])
    saveas(tmpIRIfig,[sbase '_IRI.png'])
    saveas(tmpIRIfigZoom,[sbase '_IRI_zoom.png'])
end

%% Plot lfp across space on 1 shank 
shchans = root.lfpinfo.lfpch(root.lfpinfo.lfpShank == 3);
twin = [1164.5 1165.5];
indwin = twin*root.fs_lfp ;

ripByDepthFig = figure; hold on
set(gcf,'units','normalized','position',[0.45 0.05 0.4 0.85])
for i = 1:length(shchans)
    tmpch = find(root.lfpinfo.lfpch == shchans(i));
    plot(sess.ts(root.lfp_tsb(indwin(1):indwin(2))), root.lfp(tmpch, indwin(1):indwin(2)) + i,'k')
end

if saveFlag
    saveas(ripByDepthFig,[sbase '_lfpExampleSingleShank.png'])
end

%% Get Average Ripple
nRips = size(root.ripples,1);
wdw = round(125/1000*root.fs_lfp);
ripMap = zeros(nRips, wdw*2+1);
riplf = bandpass(root.lfp(chan,:), [150, 250], root.fs_lfp);

for i = 1:nRips
    sigInds = [find(root.lfp_tsb == root.ripples(i,2)) - wdw, find(root.lfp_tsb == root.ripples(i,2)) + wdw];
    tmpSig = riplf(sigInds(1):sigInds(2));
    ripMap(i,:) = tmpSig;
end

meanRip = mean(ripMap);
sem = rmmissing(std(ripMap)/sqrt(nRips));
ciup = meanRip + sem*1.96;
cidn = meanRip - sem*1.96;
plotwdw = (-wdw:wdw) / root.fs_lfp * 1000;

avgRipFig = figure; hold on
plot(plotwdw,ripMap(1,:),'Color',[.7 .7 .7])
plot(plotwdw,ripMap(1:round(nRips/25):nRips,:),'Color',[.7 .7 .7],'HandleVisibility','off')
plot(plotwdw, meanRip, 'r')
patch([plotwdw,fliplr(plotwdw)],[cidn,fliplr(ciup)],'r','FaceAlpha',0.5,'EdgeColor','none')

xlabel('Time to ripple peak (ms)'); ylabel('150-250Hz Amplitude (mV)')

legend({'Individual Ripples','Average','95% CI'})
title(['Shank ', num2str(root.lfpinfo.lfpShank(chan)), ', Elec ' num2str(root.lfpinfo.lfpch(chan))])
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag 
    saveas(avgRipFig,[sbase '_AvgRipple.png'])
end

%% Batch plot units relative to ripples

mkdir("ripplePlots_good")
cd('ripplePlots_good')

for cc = 1:length(root.good)
    [~,~,tmpfig] = plot_frXripple(root,root.good(cc),150);
    tmpRipMod = get_RipMod(root,root.good(cc),1000);
    legend(['RippleMod ' num2str(tmpRipMod)])
    tmpind = find(root.info.cluster_id == cc);
    saveas(tmpfig, ['Unit' num2str(root.good(cc)) '_Shank' num2str(root.info.shankID(tmpind)), '_Depth', num2str(root.info.depth(tmpind)) '_ripRaster'], 'png')
    close(tmpfig)
end
cd('..')

%% Batch z-score ripple modulation and waterfall
chan = root.uPSDMax(2,1);
root.ripples = get_ripples(root,chan,sess,3,5);
% root = get_layerUnits(root,100);

nShank = numel(unique(root.info.shankID));

for i = 1:nShank
    tmpGood = find(root.goodind & root.info.shankID == i-1 & root.info.lyrID == 1);
    tmpPeriSpikes = [];
    for j = 1:length(tmpGood)
        [tmpPeriSpikes(j,:),bins] = plot_frXripple(root,root.info.cluster_id(tmpGood(j)),150,2.5,0);
    end
    figure;
    tmpMax = max(tmpPeriSpikes,[],2);
    [~,maxSort] = sort(tmpMax);
    tmpPeriSpikesSort = tmpPeriSpikes(maxSort,:);

    imagesc(tmpPeriSpikesSort)
end

%% Batch get metrics for good cells

nUnits = length(root.good);

for i = 1:nUnits
    cc = root.good(i);

    % [SI(i),uFR(i),pkFR(i)] = get_SI(root,cc,sess);
    [SI(i),pkFR(i),uFR(i)] = get_PF(root,cc,sess);
    ripMod(i) = get_RipMod(root,cc,1000);
    [~,~,vMdl] = plot_frXvel(root,cc,sess,2,0);
    vCorr(i) = vMdl.r;
end

%% Plot SI vs RipMod
ca1units = root.info.shankID(root.goodind) == 2 | root.info.shankID(root.goodind) == 3;
subunits = root.info.shankID(root.goodind) == 0 | root.info.shankID(root.goodind) == 1;
actunits = root.good;     % root.good or subunits or ca1units

mdl = fitlm(SI(actunits),ripMod(actunits));
ys = predict(mdl,SI(actunits)');
[r,p] = corrcoef(SI(actunits)',ys,'Rows','complete');
r = r(2,1);
p = p(2,1);
b = mdl.Coefficients{2,1}; 
mdlparams.r = r;
mdlparams.p = p;
mdlparams.b = b;
mdlparams.yint = predict(mdl,0);

siVrmFig = figure; hold on
plot(SI(actunits),ripMod(actunits),'k.')
plot(SI(actunits),ys,'r','LineWidth',2)
xlabel('Spatial Info. (bits/spike)'); ylabel('Ripple Modulation')
% title(['Unit ' num2str(unit)])
set(gca,'FontSize',12,'FontName','Arial')

ylims = ylim;
ylim([0 ylims(2)])
ylims = ylim;
xlims = xlim;
xlim([0 xlims(2)])
xlims = xlim;

% text(xlims(2) - .3*diff(xlims), ylims(2)-.1*diff(ylims), ['R = ' num2str(r, 3)], 'FontSize', 12)
% text(xlims(2) - .3*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(p, 3)], 'FontSize', 12)
text(xlims(2) - .3*diff(xlims), ylims(2)-.2*diff(ylims), ['slope = ' num2str(b, 3)], 'FontSize', 12)
text(xlims(2) - .3*diff(xlims), ylims(2)-.25*diff(ylims), ['y-int = ' num2str(mdlparams.yint, 3)], 'FontSize', 12)

%% Plot vCorrelation vs RipMod
mdl = fitlm(vCorr(actunits),ripMod(actunits));
ys = predict(mdl,vCorr(actunits)');
[r,p] = corrcoef(vCorr(actunits)',ys,'Rows','complete');

r = r(2,1);
p = p(2,1);
b = mdl.Coefficients{2,1};
mdlparams.r = r;
mdlparams.p = p;
mdlparams.b = b;
mdlparams.yint = predict(mdl,0);

vcVrmFig = figure; hold on
plot(vCorr(actunits),ripMod(actunits),'k.')
plot(vCorr(actunits),ys,'r','LineWidth',2)
xlabel('Spike-Velocity Tuning Corr'); ylabel('Ripple Modulation')
% title(['Unit ' num2str(unit)])
set(gca,'FontSize',12,'FontName','Arial')

ylims = ylim;
ylim([0 ylims(2)])
ylims = ylim;
xlims = xlim;
xlim([0 xlims(2)])
xlims = xlim;

% text(xlims(2) - .3*diff(xlims), ylims(2)-.1*diff(ylims), ['R = ' num2str(r, 3)], 'FontSize', 12)
% text(xlims(2) - .3*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(p, 3)], 'FontSize', 12)
text(xlims(2) - .3*diff(xlims), ylims(2)-.2*diff(ylims), ['slope = ' num2str(b, 3)], 'FontSize', 12)
text(xlims(2) - .3*diff(xlims), ylims(2)-.25*diff(ylims), ['y-int = ' num2str(mdlparams.yint, 3)], 'FontSize', 12)

if saveFlag
    saveas(siVrmFig,[root.name '_Shank' num2str(root.info.shankID(chan)) '_ripModVSI.png'])
    saveas(vcVrmFig,[root.name '_Shank' num2str(root.info.shankID(chan)) '_ripModVvCorr.png'])
end

%% Plot many LFP traces

nChans = height(root.lfpinfo);
nShank = numel(unique(root.lfpinfo.lfpShank));
twin = [1148.25 1149.25];
indwin = twin.*root.fs_lfp; % - root.lfp_tsb(1);
lfpMax = prctile(root.lfp,99,'all');
cmapcool = cool(nShank);

lfpMap = zeros(nChans, diff(indwin));
for i = 1:nChans
    lfpMap(i,:) = root.lfp(i, indwin(1):indwin(2)-1) ./ lfpMax;
end

ripXDepthXShankFig = figure; 
set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6])

for j = 1:nShank
    subplot(1,4,j)
    hold on
    ct = 0;
    shchans = find(root.lfpinfo.lfpShank == j-1);
    for i = 1:length(shchans)
        plot(sess.ts(root.lfp_tsb(indwin(1):indwin(2)-1)), lfpMap(shchans(i),:) ./ lfpMax + ct, 'Color',cmapcool(j,:))
        if shchans(i) == root.uPSDMax(2,j)
            plot(sess.ts(root.lfp_tsb(round(diff(indwin)/2+indwin(1)))), ct+0.5, '*', 'Color', cmapcool(j,:))
        end
        ct = ct + 1;
    end
    title(['Shank ' num2str(j-1)])
    yticklabels([])
    set(gca,'FontSize',12,'FontName','Arial')

end

if saveFlag; saveas(ripXDepthXShankFig,[root.name '_lfpExampleMultishank.png']); end

%% Analyze 1-2 sites/shank

chans = root.uPSDMax(2,:);
nChans = length(chans);
cmapcool = cool(nChans);

ct = 1;
catRips = [];
for i = chans
    ripStruc(ct).ripples = get_ripples(root,i,sess,4,5,[15 250]);
    catRips = [catRips; ripStruc(ct).ripples(:,2), -1+ct+zeros(size(ripStruc(ct).ripples,1),1)];
    ct = ct+1;
end

%% Plot all detected ripples across selected channels

for i = 1:nChans
    riplf(i,:) = bandpass(root.lfp(root.uPSDMax(2,i),:), [150, 250], root.fs_lfp);
end

bestChan = root.uPSDMax(2,4);
bnsz = round(50/1000*root.fs_lfp);
tmpbns = root.lfp_tsb(1):bnsz:root.lfp_tsb(end);
binRips = histcounts(catRips(:,1),tmpbns);
cmapcool = cool(nChans);

figure; hold on
set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6])

% plot(sess.ts(root.lfp_tsb), root.lfp(bestChan,:),'k')

for i = 1:length(chans)
    plot(sess.ts(root.lfp_tsb), root.lfp(root.uPSDMax(2,i),:) + (i-1), 'Color', cmapcool(i,:))
    plot(sess.ts(root.lfp_tsb), riplf(i,:) + (i*0.5)+3.5,'Color',cmapcool(i,:), 'HandleVisibility','off')
    tmprips = ripStruc(i).ripples;
    plot(sess.ts(tmprips(:,2)),6+i*0.1+zeros(1,length(tmprips(:,2))),'*','Color',cmapcool(i,:), 'HandleVisibility','off')
    plot(sess.ts(tmprips(:,1)),6+i*0.1+zeros(1,length(tmprips(:,1))),'|','Color',cmapcool(i,:), 'HandleVisibility','off')
    plot(sess.ts(tmprips(:,3)),6+i*0.1+zeros(1,length(tmprips(:,3))),'|','Color',cmapcool(i,:), 'HandleVisibility','off')
end
% plot(sess.ts(tmpbns(1:end-1)),binRips/max(binRips)+2,'k')
plot(sess.ts,sess.velshft/max(sess.velshft)-2,'k')
xlabel('Time (sec)')
legCell = {'Sh0','Sh1','Sh2','Sh3','Velocity'};
legend(legCell)

%% Compare histogram of ripples between two sites/shanks
compParam = 2;  % 1 = rip start, 2 = rip peak, 3 = rip end

nearStruc = struct;

for i = 1:nChans
    for j = 1:nChans
        nearStruc(i).nearInds(j,:) = compare_ripples_times(ripStruc(i).ripples(:,compParam), ripStruc(j).ripples(:,compParam), root.fs_lfp, 0.5);
    end
end
nearStruc4 = nearStruc;

bnsz = 0.001;
wnlen = 0.05;
bins = -wnlen:bnsz:wnlen;

rippleHistogramFig = figure; ct = 1;
set(gcf,'units','normalized','position',[0.01 0.01 0.95 0.9])
for i = 1:nChans
    for j = 1:nChans
        if i == j
            ct = ct + 1;
            continue
        else
            subplot(nChans,nChans,ct)
            hold on
            nearRips = nearStruc(i).nearInds(j,:);
            matchRips = find(~isnan(nearRips));
            rDist = sess.ts(nearRips(matchRips)) - sess.ts(ripStruc(i).ripples(matchRips,compParam));
            binnedRipCounts = histcounts(rDist,bins);
            plot([0 0], [0 max(binnedRipCounts)+10],'k--','HandleVisibility','off')
            bar((bins(1:end-1)+0.5*bnsz),binnedRipCounts,'b');
            title(['Sh' num2str(i-1) ' minus Sh' num2str(j-1)]);
            % legend([num2str(length(matchRips)) ' Matches'])
            ct = ct + 1;
        end
    end
end
han=axes(rippleHistogramFig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Time to ripple peak (sec)')
ylabel(han,'Count')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag; saveas(rippleHistogramFig,[root.name '_rippleTimingHisto.png']); end

%% Layer histograms atop each other

layerHistoFig = figure;
set(gcf,'units','normalized','position',[0.01 0.3 0.95 0.3])
for i = 1:nChans
    subplot(1,nChans,i); hold on
    plot([0 0], [0 0.5],'k--','HandleVisibility','off')
    legCell = {};
    maxVal = 0;
    for j = 1:nChans
        if i == j
            continue
        else
            legCell = [legCell ['Sh' num2str(j-1)]];
            nearRips = nearStruc(i).nearInds(j,:);
            matchRips = find(~isnan(nearRips));
            rDist = sess.ts(nearRips(matchRips)) - sess.ts(ripStruc(i).ripples(matchRips,compParam));
            plot(0,0,'Color',cmapcool(j,:))
            h = histogram(rDist,bins,'EdgeColor',cmapcool(j,:),'DisplayStyle','stairs','Normalization','probability','HandleVisibility','off');
            if maxVal < max(h.Values); maxVal = max(h.Values); end
        end
    end
    xlabel(['Time (ms) to Peak on Sh' num2str(i-1)]);
    ylim([0 maxVal])
    legend(legCell)
end
han=axes(layerHistoFig,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Normalized Count Probability')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag; saveas(layerHistoFig,[root.name '_rippleTimingHistoCollapse.png']); end

%% Cross correlate ripple signals
compParam = 2;  % 1 = rip start, 2 = rip peak, 3 = rip end
winlen = 0.05;

rippleXCorrFig = figure; ct = 1;
set(gcf,'units','normalized','position',[0.01 0.01 0.95 0.9])
for i = 1:nChans
    for j = 1:nChans
        if i == j
            ct = ct + 1;
            continue
        else
            subplot(nChans,nChans,ct)
            hold on
            rip1vec = histcounts(ripStruc(i).ripples(:,compParam),root.lfp_tsb);
            rip2vec = histcounts(ripStruc(j).ripples(:,compParam),root.lfp_tsb);
            [ripcorr,lags] = xcorr(rip2vec, rip1vec, round(root.fs_lfp*winlen), 'normalized');
            plot([0 0], [0 max(ripcorr)],'k--')
            plot(lags/root.fs_lfp,ripcorr,'b');
            title(['Sh' num2str(i-1) ' minus Sh' num2str(j-1)]);
            xlim([-winlen winlen])
            ct = ct + 1;
        end
    end
end
han=axes(rippleXCorrFig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Time to ripple peak (sec)')
ylabel(han,'Normalized AutoCorr')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag; saveas(rippleXCorrFig,[root.name '_rippleTimingXCorr.png']); end

%% Layer histograms atop each other

layerXCorrFig = figure;
set(gcf,'units','normalized','position',[0.01 0.3 0.95 0.3])
for i = 1:nChans
    subplot(1,nChans,i); hold on
    plot([0 0], [0 1],'k--','HandleVisibility','off')
    legCell = {};
    maxCorr = 0;
    for j = 1:nChans
        if i == j
            continue
        else
            legCell = [legCell ['Sh' num2str(j-1)]];
            rip1vec = histcounts(ripStruc(i).ripples(:,compParam),root.lfp_tsb);
            rip2vec = histcounts(ripStruc(j).ripples(:,compParam),root.lfp_tsb);
            [ripcorr,lags] = xcorr(rip2vec, rip1vec, round(root.fs_lfp*winlen), 'normalized');
            plot(lags/root.fs_lfp*1000,ripcorr,'Color',cmapcool(j,:));
            if maxCorr < max(ripcorr); maxCorr = max(ripcorr); end
        end
    end
    xlabel(['Time (ms) to Peak on Sh' num2str(i-1)]);
    ylim([0 maxCorr])
    legend(legCell)
end
han=axes(layerXCorrFig,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Normalized Count Probability')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag; saveas(layerXCorrFig,[root.name '_rippleTimingXCorrCollapse.png']); end

%% Concatenate ripple timing comparisons across sessions

% cd('D:\Data\Kelton\analyses\KW022\KW022_12162024_rec_D1_RLat1')
% load('KW022_12162024_session.mat')
% load('KW022_12162024_rec_D1_RLat1_root.mat')
% shCA1 = 4;
% shSub = 3;
cd('D:\Data\Kelton\analyses\KW022\KW022_12172024_rec_D2_RMed1')
load('KW022_12172024_session.mat')
load('KW022_12172024_rec_D2_RMed1_root.mat')
shCA1 = 4;
shSub = 3;
% cd('D:\Data\Kelton\analyses\KW022\KW022_12202024_rec_D4_LMed1')
% load('KW022_12202024_session.mat')
% load('KW022_12202024_rec_D4_LMed1_root.mat')
% shCA1 = 4;
% shSub = 3;

shCA1Ripples = get_ripples(root,root.uPSDMax(2,shCA1),sess,4,5,[15 250]);
shSuBRipples = get_ripples(root,root.uPSDMax(2,shSub),sess,2,2,[15 250]);

bnsz = 0.001;
wnlen = 0.05;
bins = -wnlen:bnsz:wnlen;

% rDistT1 = getRDistT(root,sess,shCA1Ripples(:,2),shSuBRipples(:,2));
% rDistT2 = getRDistT(root,sess,shCA1Ripples(:,2),shSuBRipples(:,2));
% rDistT4 = getRDistT(root,sess,shCA1Ripples(:,2),shSuBRipples(:,2));
rDistTAll = [rDistT1, rDistT2, rDistT4];

figure; hold on
binnedRipCounts = histcounts(rDistTAll,bins);
plot([0 0], [0 max(binnedRipCounts)+10],'k--','HandleVisibility','off')
bar((bins(1:end-1)+0.5*bnsz),binnedRipCounts,'b');
title('CA1 Vs CA1');
xlabel('Time to CA1 Prox ripple peak (s)')
ylabel('Matched Ripple Count')
legend('Sub Thresh = [2,2]')

%% Functions

function [rDistT] = getRDistT(root,sess,rInds1,rInds2)

nearInds = compare_ripples_times(rInds1, rInds2, root.fs_lfp, 0.5);
matchRips = find(~isnan(nearInds));
rDistT = sess.ts(nearInds(matchRips)) - sess.ts(rInds1(matchRips));

end