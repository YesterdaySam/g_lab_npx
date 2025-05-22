%% Prepare data table and run parameters
parentDir = "D:\Data\Kelton\analyses\group_analyses\Subiculum_Speed"; 
subDatT = import_xldat(parentDir,"dat_include.xlsx");
cd(parentDir) 

cleanInds = subDatT.include == 0;
subDatT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;

winlen = 2;

%%

subSpeedTimes = [];
subSpeedSlope = [];
unitTypes = [];
mouse = [];
recID = [];

for i = 1:7 %height(subDatT)
    cd(subDatT.fpath{i})

    rootfile = dir("*_root.mat");
    load(rootfile.name)
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    disp(root.name)

    nShank = numel(unique(root.info.shankID));

    % root = get_layerUnits(root,100);
    try
        root.info.uType(1);
    catch
        % root = get_estCellType(root,15,5,100,1);    %Screen putative interneurons based on width and FWHM thresholds, not FR
        root = get_umapCellType(root,0);    %Screen putative interneurons based on UMAP and k-means clustering
    end
    root = get_FRVar(root,sess,0);  %Screen for overly high variance of firing rate

    for j = 1:nShank
        disp(['Shank' num2str(j-1)])
        if subDatT{i,6+j}{1} == 'sub'
            lyrUnits = find(root.info.shankID == j-1 & root.info.lyrID == 1 & root.goodind & root.info.frVar < 1.4);
            lyrPyrs  = find(root.info.shankID == j-1 & root.info.lyrID == 1 & root.goodind & root.info.uType == 1 & root.info.frVar < 1.4);
            lyrINs   = find(root.info.shankID == j-1 & root.info.lyrID == 1 & root.goodind & root.info.uType == 0 & root.info.frVar < 1.4);
            disp(['N = ' num2str(numel(lyrPyrs)) ' good pyrs and N = ' num2str(numel(lyrINs)) ' good INs'])
            disp(['N = ' num2str(sum(root.info.frVar(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1) < 1.4)),...
                ' stable units and N = ' num2str(sum(root.info.frVar(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1)>=1.4)) ' high variance units'])

            for k = 1:length(lyrUnits)
                cc = root.info.cluster_id(lyrUnits(k));
                [~,tmpTime,tmpCorr] = shiftVelCoding(root,cc,sess,winlen,2,'pearson',0); % window +/-0.5sec, 2cm/s binsize, corrtype, no plot
                subSpeedTimes = [subSpeedTimes; tmpTime];
                subSpeedSlope = [subSpeedSlope; tmpCorr > 0];
            end
            unitTypes = [unitTypes; root.info.uType(lyrUnits)];
            mouse = [mouse; repmat(str2num(subDatT.mouse{i}(end-1:end)), length(lyrPyrs),1)];
            recID = [recID; ones(length(lyrPyrs),1)*i];
        end
    end
end
% subSpeedTimes = subSpeedTimes * 1000;

%%
winlen = 2;
% unitTypes = logical(unitTypes);
% subSpeedSlope = logical(subSpeedSlope);
% 
tShiftBnsz = 0.050;
tShifts = -winlen:tShiftBnsz:winlen;
% 
binTPyrs = histcounts(subSpeedTimes(unitTypes),tShifts);
binTPyrs = binTPyrs / sum(binTPyrs);

figure;
subplot(2,1,1)
bar(1000*tShifts(1:end-1)+tShiftBnsz*0.5,binTPyrs,'b')
xlabel('Optimum Coding lag (ms)'); ylabel('Counts (Norm)')
title('Pyramidals')
legend(['Mean: ' num2str(1000*mean(subSpeedTimes(unitTypes)),3)])
set(gca,'FontSize',12,'FontName','Arial')

binTINs = histcounts(subSpeedTimes(~unitTypes),tShifts);
binTINs = binTINs / sum(binTINs);

subplot(2,1,2)
bar(1000*tShifts(1:end-1)+tShiftBnsz*0.5,binTINs,'r')
xlabel('Optimum Coding lag (ms)'); ylabel('Counts (Norm)')
legend(['Mean: ' num2str(1000*mean(subSpeedTimes(~unitTypes)),3)])
title('Interneurons')
set(gca,'FontSize',12,'FontName','Arial')

%%
binTPosPyrs = subSpeedTimes(unitTypes & subSpeedSlope);
% binTPosPyrsN = binTPosPyrs / sum(binTPosPyrs);
binTNegPyrs = subSpeedTimes(unitTypes & ~subSpeedSlope);
% binTNegPyrsN = binTNegPyrs / sum(binTNegPyrs);

figure; hold on
plot([mean(binTPosPyrs) mean(binTPosPyrs)],[0 0.1],'b--')
plot([mean(binTNegPyrs) mean(binTNegPyrs)],[0 0.1],'r--')
histogram(binTPosPyrs,tShifts,'EdgeColor','b','DisplayStyle','stairs','LineWidth',2,'Normalization','probability')
histogram(binTNegPyrs,tShifts,'EdgeColor','r','DisplayStyle','stairs','LineWidth',2,'Normalization','probability')
xlabel('Optimum Coding lag (s)'); ylabel('Counts (Norm)')
legend({['Mean lag (+FR): ' num2str(mean(binTPosPyrs))],['Mean lag (-FR): ' num2str(mean(binTNegPyrs))]},'location','northwest')
title('Pyramidals')
set(gca,'FontSize',12,'FontName','Arial')

%% Plot individual examples
unit = 188;

bnsz = 0.5;
tbins = sess.ts(1):bnsz:sess.ts(end);
spkinds = root.tsb(root.cl == unit);
binfr = histcounts(sess.ts(spkinds),tbins)/bnsz;
binfr = smooth(binfr);
binfr = binfr./max(binfr);
% binvel = downsample(sess.velshft,250);
binvel = smooth(sess.velshft./max(sess.velshft));

figure; hold on
plot(sess.ts,binvel,'Color',[.5 .5 .5],'LineWidth',1)
plot(tbins(1:end-1)+bnsz/2,binfr,'r')
xlabel('Time (s)'); ylabel('Normalized Value')
legend({'Velocity','Smoothed FR'})

%% Cross Correlation method
winlen = 1*sess.samprate;

% [ccorr,clags] = xcorr(binvel(2:end),binfr);
[ccorr,clags] = xcorr(sess.velshft,histcounts(sess.ts(spkinds),sess.ts),winlen);
maxlag = clags(ccorr == max(ccorr));
figure; hold on
plot(clags/sess.samprate*1000,ccorr)
plot([0 0],[0 max(ccorr) + 1000],'k--')
plot(maxlag/sess.samprate*1000,max(ccorr),'r*')
xlim([-winlen/sess.samprate*1000 winlen/sess.samprate*1000])
xlabel('Lag (ms)'); ylabel('Cross Correlation')


%% Subiculum Theta Phase code test
ps = [];
angs = [];
mrls = [];

% Pre-calculate theta at s.pyr. center, based on 150-250Hz PSD peak across shanks
for i = 1:4
    thetaLFPSh = bandpass(root.lfp(root.uPSDMax(2,i),:),root.bands(1,:), root.fs_lfp);
    root.thEnv(i,:) = hilbert(thetaLFPSh);  % Operates along rows
end

for i = 1:length(root.good)
    cc = root.good(i);
    tmpsh = root.info.shankID(find(root.info.cluster_id == cc));
    tmp_stats = plot_thetaMod(root,cc,tmpsh+1,2*pi/18,0);    %Using default 20 deg bins
    ps(i) = tmp_stats.p;
    angs(i) = tmp_stats.ang;
    mrls(i) = tmp_stats.mrl;
end

figure; histogram(rad2deg(angs(ps < 0.05)),rad2deg(-pi:2*pi/18:pi))
hold on; histogram(rad2deg(angs(ps > 0.05)),rad2deg(-pi:2*pi/18:pi),'EdgeColor','r','DisplayStyle','stairs')
legend({'p < 0.05','p > 0.05'})
xlabel('Theta Phase (0 Trough)'); ylabel('Count')
set(gca,'FontSize',12,'FontName','Arial')

figure; histogram(mrls(ps < 0.05),0:max(mrls)/20:max(mrls))
hold on; histogram(mrls(ps > 0.05),0:max(mrls)/20:max(mrls),'EdgeColor','r','DisplayStyle','stairs')
legend({'p < 0.05','p > 0.05'})
xlabel('Mean Resultant Length'); ylabel('Count')
set(gca,'FontSize',12,'FontName','Arial')

%% Example Units
cc = 217;
tmpsh = root.info.shankID(find(root.info.cluster_id == cc));
figure; plot(sess.ts(root.lfp_tsb),bandpass(root.lfp(root.uPSDMax(2,tmpsh+1),:),root.bands(1,:),root.fs_lfp),'r')
hold on
spk_tsb = root.tsb(root.cl == cc);
plot(sess.ts(spk_tsb),ones(length(spk_tsb),1)-0.8,'k|')

tmp_stats = plot_thetaMod(root,cc,tmpsh+1);