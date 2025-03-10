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

saveFlag = 0;

% chan = find(root.lfpinfo.lfpch == 336);
chan = root.uPSDMax(2,4);

root.ripples = get_ripples(root,chan,sess);

rawlf = root.lfp(chan,:);
riplf = bandpass(rawlf, [150, 250], root.fs_lfp);

%% Optional plotting

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
    saveas(ripByDepthFig,[sbase '_RippleByDepthExample.png'])
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
    [~,~,tmpfig] = plot_frXripple(root,root.good(cc));
    tmpRipMod = get_RipMod(root,root.good(cc),1000);
    legend(['RippleMod ' num2str(tmpRipMod)])
    tmpind = find(root.info.cluster_id == cc);
    saveas(tmpfig, ['Unit' num2str(root.good(cc)) '_Shank' num2str(root.info.shankID(tmpind)), '_Depth', num2str(root.info.depth(tmpind)) '_ripRaster'], 'png')
    close(tmpfig)
end
cd('..')

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
secs = [635 655];
inds = secs.*root.fs_lfp - root.lfp_tsb(1);
lfpMax = prctile(root.lfp,99,'all');

% figure; hold on
% ct = 0;
% for i = 1:6
%     plot(sess.ts(root.lfp_tsb(inds(1):inds(2))), root.lfp(i, inds(1):inds(2)) ./ lfpMax + ct, 'k')
%     ct = ct + 1;
% end
% for i = 13:18
%     plot(sess.ts(root.lfp_tsb(inds(1):inds(2))), root.lfp(i, inds(1):inds(2)) ./ lfpMax + ct, 'k')
%     ct = ct + 1;
% end

lfpMap = zeros(nChans, diff(inds));
for i = 1:nChans
    lfpMap(i,:) = root.lfp(i, inds(1):inds(2)-1) ./ lfpMax;
end

[~,sortInds] = sort(root.lfpinfo.lfpShank);
lfpMap = lfpMap(sortInds,:);

figure; 
set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6])

subplot(1,4,1)
ct = 0;
hold on
for i = 1:12
    plot(sess.ts(root.lfp_tsb(inds(1):inds(2)-1)), lfpMap(i,:) ./ lfpMax + ct, 'k')
    ct = ct + 1;
end
title('Shank 0')
yticklabels([])
set(gca,'FontSize',12,'FontName','Arial')

subplot(1,4,2)
ct = 0;
hold on
for i = 13:24
    plot(sess.ts(root.lfp_tsb(inds(1):inds(2)-1)), lfpMap(i,:) ./ lfpMax + ct, 'k')
    ct = ct + 1;
end
title('Shank 1')
yticklabels([])
set(gca,'FontSize',12,'FontName','Arial')

subplot(1,4,3)
ct = 0;
hold on
for i = 25:36
    plot(sess.ts(root.lfp_tsb(inds(1):inds(2)-1)), lfpMap(i,:) ./ lfpMax + ct, 'k')
    ct = ct + 1;
end
title('Shank 2')
yticklabels([])
set(gca,'FontSize',12,'FontName','Arial')

subplot(1,4,4)
ct = 0;
hold on
for i = 37:nChans
    plot(sess.ts(root.lfp_tsb(inds(1):inds(2)-1)), lfpMap(i,:) ./ lfpMax + ct, 'k')
    ct = ct + 1;
end
title('Shank 3')
yticklabels([])
set(gca,'FontSize',12,'FontName','Arial')

saveas(gcf,[root.name '_lfpExample.png'])

%% Analyze 1-2 sites/shank

chans = root.uPSDMax(2,:);
nChans = length(chans);
ct = 1;
catRips = [];
for i = chans
    ripStruc(ct).ripples = get_ripples(root,i,sess);
    catRips = [catRips; ripStruc(ct).ripples(:,2), -1+ct+zeros(size(ripStruc(ct).ripples,1),1)];
    ct = ct+1;
end

%% Plot all detected ripples
bestChan = 1;
bnsz = round(50/1000*root.fs_lfp);
tmpbns = root.lfp_tsb(1):bnsz:root.lfp_tsb(end);
binRips = histcounts(catRips(:,1),tmpbns);

figure; hold on
set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6])

plot(sess.ts(root.lfp_tsb), root.lfp(chans(bestChan),:),'k')

for i = 1:length(chans)
    tmprips = ripStruc(i).ripples(:,2);
    plot(sess.ts(tmprips),3+i*0.1+zeros(1,length(tmprips)),'*')
end
plot(sess.ts(tmpbns(1:end-1)),binRips/max(binRips)+2)

%%
bestChan = 1;
nRips = size(ripStruc(bestChan).ripples,1);
nearRips = zeros(nRips,nChans);
for i = 1:nRips
    pkInd = ripStruc(bestChan).ripples(i,2);
    % nearRips(i,1) = pkInd;
    for j = 1:nChans
        [tmpMin, tmpMinInd] = min(abs(pkInd - ripStruc(j).ripples(:,2)));
        if tmpMin < 1*root.fs_lfp     % 1 second
            nearRips(i,j) = ripStruc(j).ripples(tmpMinInd,2);
        else
            nearRips(i,j) = NaN;
        end
    end
end

%%
rDist = nearRips(:,2:nChans) - nearRips(:,1);
figure; hold on
histogram(rDist)

bins = -root.fs_lfp:2.5:root.fs_lfp;
binnedRipCounts = histcounts(reshape(rDist,1,[]),bins);

figure; hold on
plot(rDist(:,1)/root.fs_lfp,1:nRips,'k|')
for i = 3:2:nChans  
    plot(rDist(:,i)/root.fs_lfp,1:nRips,'|')
end
plot((bins(1:end-1)+1.25)/root.fs_lfp,binnedRipCounts,'r')
xlim([-0.02 0.02])
xlabel(['Time to reference ripple chan ' num2str(chans(bestChan))])
ylabel('Ripple #')