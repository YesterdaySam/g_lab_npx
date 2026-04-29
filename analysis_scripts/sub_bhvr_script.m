%% sub_bhvr_script

parentDir = 'D:\Data\Kelton\analyses\AT_behavior';

overwriteFlag = 1;

cd(parentDir)
dirlist = dir;
dirlist = dirlist([dirlist.isdir]);
for i = 1:numel(dirlist)
    useDir(i) = contains(dirlist(i).name,'AT0');
end
dirlist = dirlist(useDir);   % Remove parent dirs

%% Import behavior

for i = 1:length(dirlist)
    cd(fullfile(dirlist(i).folder,dirlist(i).name));

    sessmade = dir('*session.mat');
    if ~isempty(sessmade) & overwriteFlag == 0
        disp(['Not a valid directory or session already created for ' dirlist(i).name])
        cd(parentDir)
        continue
    end

    try
        sess = importBhvr(fullfile(dirlist(i).folder,dirlist(i).name));
        sess.maxPos = 2;    % Account for AT23-25 longer belts
        sess = getErrorTrials(sess);
    catch
        cd(parentDir)
        continue
    end

    try
        [~,sess.velXpos] = plot_trialvel(sess,0.01,0);
        [~,sess.velXrwd] = get_velXrwd(sess,0.25,5,0);
    catch
        disp('No bnvel or velxrwd assigned to sess file')
    end

    sbase = sess.name(1:12);
    save([sbase, '_session'], 'sess','-v7.3')
    disp(['Session created and saved for ' sbase])

    if isempty(sess.valTrials)
        continue
    end

    [tmpedges1, tmpbnvel, fig_trialvel,fig_velavg] = plot_trialvel(sess);  % default 0.01m binsize
    saveas(fig_trialvel,[sbase, '_trialvelocity'],'png')
    % saveas(fig_velavg,[sbase, '_avgvelocity'],'png')

    [fig_lickraster, fig_lickavg] = plot_lickraster_rwd(sess, 5, 0.1);
    saveas(fig_lickraster,[sbase, '_lickraster_rwdalign'],'png')
    saveas(fig_lickavg,[sbase, '_lickaverage'],'png')

    [tmpedges2, ~, tmpbnlck, fig_lickpos, fig_licktrialavg] = plot_lickpos(sess);   % Don't run this for D1 - random acclimation (# laps ~= # rewards)
    saveas(fig_lickpos,[sbase, '_lickraster'],'png')
    % saveas(fig_licktrialavg,[sbase, '_lickpos_average'],'png')

    fig_vel_lick = plot_vel_lck(sess, tmpbnvel, 0.01, 0.03);
    saveas(fig_vel_lick,[sbase, '_vel_lck'],'png')

    close all

    plotBhvrCompact(sess,fullfile(dirlist(i).folder,dirlist(i).name),1)

    disp(['Finished basic import and visualization of ', sbase])

end

disp(['Finished basic import and visualization of all files for ' parentDir])

%% Fix early reset sessions
% sessOrig = sess;

lapLengths = sess.pos(sess.lapend) - sess.pos(sess.lapstt);
shortLaps = lapLengths < 1.9;
addLap = zeros(size(lapLengths));

for i = 2:length(lapLengths)-1
    if ~shortLaps(i)
        continue
    end
    tmpLen = sess.pos(sess.lapend(i));
    preLen = sess.pos(sess.lapend(i-1));
    
    if tmpLen + preLen > 2.1
        continue
    end
    addLap(i) = true;
    sess.pos(sess.lapstt(i):sess.lapend(i)) = sess.pos(sess.lapstt(i):sess.lapend(i)) + sess.pos(sess.lapend(i-1));
end


addLap              = logical(addLap);
sess.lapstt(addLap) = [];
sess.lapend         = [sess.lapstt(2:end) - 1; sess.ind(end)];     %Use last ts as last lap end
sess                = getErrorTrials(sess);
sess.nlaps          = length(sess.lapstt);

saveSess(sess)

%% Set parameters

r1pos = 0.4;    % 20-40 cm
r2pos = 1.3;      % 0-20 cm
zLen  = 30;
vspbnsz = 0.01;
sttendSz = 20;  % N trials for averaging early vs late metrics in session
mExp = [5, 12, 14, 21, 24, 25];
mCtl = [8 9 11 13 15 22 23];

saveFlag = 0;
sdir = 'D:\Data\Kelton\analyses\AT_behavior\groupAnalyses\az_30cm_lap100';
sName = 'AT_Dreadds_';

%% Derive Behavioral scores for each session

clear mID rDay gExp nLap bhvr lsiSttEnd lapScoreSttEnd vCorSttEnd
for i = 1:length(dirlist)
    cd(fullfile(dirlist(i).folder,dirlist(i).name));
    
    sbase = dirlist(i).name;
    mID(i) = str2num(dirlist(i).name(3:5));
    rDay(i) = str2num(dirlist(i).name(end));
    gExp(i) = ~isempty(find(mExp == mID(i),1));

    sessF = dir("*_session.mat");
    load(sessF.name)
    nLap(i) = length(sess.valTrials);
    lrange = [1:sttendSz; min(100, nLap(i))-sttendSz+1:min(100, nLap(i))];

    % LSI
    [prepLckHz, bhvr(i).lapLSI] = get_lickDiscrim(sess,[r1pos r2pos]*100,zLen);
    lsiFig = plot_lickDiscrim(sess,[r1pos r2pos]*100,zLen,10);
    lsiSttEnd(i,:) = [mean(bhvr(i).lapLSI(lrange(1,:))) mean(bhvr(i).lapLSI(lrange(2,:)))];

    % Operant success
    bhvr(i).lapScore = get_trialSuccess(sess);
    opScoreFig = plot_trialScore(sess);
    lapScoreSttEnd(i,:) = [mean(bhvr(i).lapScore(lrange(1,:))) mean(bhvr(i).lapScore(lrange(2,:)))];

    % Velocity
    bhvr(i).vCorXLap = corr(sess.velXpos', mean(sess.velXpos(end-19:end,:))','rows','complete');
    tmpVCorrFig = figure; hold on; plot(sess.valTrials,bhvr(i).vCorXLap,'k-o',sess.valTrials,smooth(bhvr(i).vCorXLap,5),'r')
    ylim([-0.5 1]); ylabel('Velocity Corr. to last 20 laps'); xlabel('Lap #'); 
    set(gca,'FontSize',16,'FontName','Arial')
    vCorSttEnd(i,:) = [mean(bhvr(i).vCorXLap(lrange(1,:))) mean(bhvr(i).vCorXLap(lrange(2,:)))];

    if saveFlag
        fsave(lsiFig,[sbase '_lsiXtrial'],1,0,0)
        fsave(opScoreFig,[sbase '_operantScoreXtrial'],1,0,0)
        fsave(tmpVCorrFig,[sbase '_velocityCorrxtrial'],1,0,0)
    end
    close all
end

cd(sdir)

if saveFlag
    save('AT_bhvr_vars_30cm_lap100','mID','rDay','gExp','nLap','bhvr','lsiSttEnd','lapScoreSttEnd','vCorSttEnd');
end

%% Compare delta First 20 vs Last 20 laps per condition

expMice = gExp == 1;
nMice = [numel(unique(mID))-numel(mExp), numel(mExp)];  % [ctl, exp]
nDays = 4; % Ignoring shift day for now
vColors = [.35 .35 .35; 1 .25 .25];

lsiDlt = lsiSttEnd(:,2) - lsiSttEnd(:,1);
lpsDlt = lapScoreSttEnd(:,2) - lapScoreSttEnd(:,1);
vlcDlt = vCorSttEnd(:,2) - vCorSttEnd(:,1);

for i = 1:nDays
    useDays = rDay == i;
    ctlLSI_frst(:,i) = lsiSttEnd(useDays & ~expMice,1);
    ctlLSI_last(:,i) = lsiSttEnd(useDays & ~expMice,2);
    expLSI_frst(:,i) = lsiSttEnd(useDays & expMice,1);
    expLSI_last(:,i) = lsiSttEnd(useDays & expMice,2);
    ctlLpS_frst(:,i) = lapScoreSttEnd(useDays & ~expMice,1);
    ctlLpS_last(:,i) = lapScoreSttEnd(useDays & ~expMice,2);
    expLpS_frst(:,i) = lapScoreSttEnd(useDays & expMice,1);
    expLpS_last(:,i) = lapScoreSttEnd(useDays & expMice,2);
    ctlVlC_frst(:,i) = vCorSttEnd(useDays & ~expMice,1);
    ctlVlC_last(:,i) = vCorSttEnd(useDays & ~expMice,2);
    expVlC_frst(:,i) = vCorSttEnd(useDays & expMice,1);
    expVlC_last(:,i) = vCorSttEnd(useDays & expMice,2);

    [~,ps(i).lsi_delta_expVctl] = ttest2(lsiDlt(useDays & expMice), lsiDlt(useDays & ~expMice));
    [~,ps(i).lps_delta_expVctl] = ttest2(lpsDlt(useDays & expMice), lpsDlt(useDays & ~expMice));
    [~,ps(i).vlc_delta_expVctl] = ttest2(vlcDlt(useDays & expMice), vlcDlt(useDays & ~expMice));
    [~,ps(i).lap_comps_expVctl] = ttest2(nLap(useDays & expMice), nLap(useDays & ~expMice));

    ctlDays(:,i) = rDay(useDays & ~expMice);
    expDays(:,i) = rDay(useDays & expMice);

    tmpLSIF = plot_barXgroup(lsiDlt(useDays & ~expMice), lsiDlt(useDays & expMice), vColors);
    ylim([-0.2 0.75]); text2bar(tmpLSIF,'\Delta LSI',ps(i).lsi_delta_expVctl);

    tmpLpSF = plot_barXgroup(lpsDlt(useDays & ~expMice), lpsDlt(useDays & expMice), vColors);
    ylim([-0.2 0.75]); text2bar(tmpLpSF,'\Delta Lap Score',ps(i).lps_delta_expVctl);

    tmpVlCF = plot_barXgroup(vlcDlt(useDays & ~expMice), vlcDlt(useDays & expMice), vColors);
    ylim([-0.2 0.75]); text2bar(tmpVlCF,'\Delta Velocity Corr',ps(i).vlc_delta_expVctl);

    tmpLapF = plot_barXgroup(nLap(useDays & ~expMice), nLap(useDays & expMice), vColors);
    ylim([0 250]); text2bar(tmpLapF,'# Laps',ps(i).lap_comps_expVctl);

    if saveFlag
        fsave(tmpLSIF,[sName '_LSI_delta_D' num2str(i)],1,0,0)
        fsave(tmpLpSF,[sName '_LpS_delta_D' num2str(i)],1,0,0)
        fsave(tmpVlCF,[sName '_VlC_delta_D' num2str(i)],1,0,0)
        fsave(tmpLapF,[sName '_nLaps_D' num2str(i)],1,0,0)
        close all
    end
end

%% Graph raw data split first/last 20 across days

lsiXtime_frstF = plot_bhvrXdays(ctlLSI_frst,expLSI_frst,ctlDays,expDays,vColors);
ylim([-0.2 1]); ylabel('LSI')
lsiXtime_lastF = plot_bhvrXdays(ctlLSI_last,expLSI_last,ctlDays,expDays,vColors);
ylim([-0.2 1]); ylabel('LSI')

lpSXtime_frstF = plot_bhvrXdays(ctlLpS_frst,expLpS_frst,ctlDays,expDays,vColors);
ylim([-0.2 1]); ylabel('Operant Success')
lpSXtime_lastF = plot_bhvrXdays(ctlLpS_last,expLpS_last,ctlDays,expDays,vColors);
ylim([-0.2 1]); ylabel('Operant Success')

vlcXtime_frstF = plot_bhvrXdays(ctlVlC_frst,expVlC_frst,ctlDays,expDays,vColors);
ylim([-0.1 1]); ylabel('Velocity Correlation')
vlcXtime_lastF = plot_bhvrXdays(ctlVlC_last,expVlC_last,ctlDays,expDays,vColors);
ylim([-0.1 1]); ylabel('Velocity Correlation')

if saveFlag
    fsave(lsiXtime_frstF,[sName 'LSI_first20_Xdays'],1,0,0)
    fsave(lsiXtime_lastF,[sName 'LSI_last20_Xdays'],1,0,0)
    fsave(lpSXtime_frstF,[sName 'LpS_first20_Xdays'],1,0,0)
    fsave(lpSXtime_lastF,[sName 'LpS_last20_Xdays'],1,0,0)
    fsave(vlcXtime_frstF,[sName 'VlC_first20_Xdays'],1,0,0)
    fsave(vlcXtime_lastF,[sName 'VlC_last20_Xdays'],1,0,0)
    % close all
end

%% Plot 1 day LSI over time

plotDay = 2;
lapMexp = nan(numel(mExp),max(nLap(rDay == plotDay & gExp)));
lapMctl = nan(numel(mCtl),max(nLap(rDay == plotDay & ~gExp)));
cte = 1;
ctc = 1;

for i = 1:length(dirlist)
    if rDay(i) ~= plotDay
        continue
    end
    if gExp(i) == 1
        vcol = 'r';
        lapMexp(cte,1:nLap(i)) = bhvr(i).lapLSI;
        cte = cte + 1;
    else
        vcol = 'k';
        lapMctl(ctc,1:nLap(i)) = bhvr(i).lapLSI;
        ctc = ctc + 1;
    end

    % plot(bhvr(i).lapLSI,vcol)

end

lapXlapsF = figure; hold on
[expUp,expDn] = get_CI(lapMexp);
[ctlUp,ctlDn] = get_CI(lapMctl);
plot_CIs(1:length(expUp),expUp,expDn,'r')
plot_CIs(1:length(ctlUp),ctlUp,ctlDn,'k')

plot(mean(lapMexp,'omitnan'),'r','LineWidth',2)
plot(mean(lapMctl,'omitnan'),'k','LineWidth',2)
ylabel('LSI'); ylim([-1 1])
xlabel('Lap #'); xlim([0.5 max(nLap(rDay == plotDay))])
title(['Day ' num2str(plotDay)])
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    fsave(lapXlapsF,[sName 'LSI_Xlaps_D' num2str(plotDay)],1,0,0)
end

%% Analyze D5 RZ shift

sdir = 'D:\Data\Kelton\analyses\AT_behavior\groupAnalyses\d5_30cm_lsi';
% clear mID rDay gExp nLap
useDays = rDay == 5;
useSess = find(useDays);
famR1 = 0.4; famR2 = 0.1;
novR1 = 1.3; novR2 = 1;
useLaps = 1:30;
ct = 1;

for i = useSess
    cd(fullfile(dirlist(i).folder,dirlist(i).name));
    
    sbase = dirlist(i).name;
    % famDat(ct).mID = str2num(dirlist(i).name(3:5));
    % famDat(ct).rDay = str2num(dirlist(i).name(end));
    % famDat(ct).gExp = ~isempty(find(mExp == mID(i),1));

    sessF = dir("*_session.mat");
    load(sessF.name)
    % famDat(ct).nLap = length(sess.valTrials);
    % lrange = [1:sttendSz; min(100, famDat(ct).nLap)-sttendSz+1:min(100, famDat(ct).nLap)];

    rwdShift = find(diff(sess.pos(sess.rwdind)) > 0.4,1);   % Find lap of reward shift
    if isfield(sess,'rwdTrials')
        rwdShift = sess.rwdTrials(rwdShift);
    end
    famInds         = [sess.ind(1) sess.lapend(rwdShift)];
    novInds         = [sess.lapstt(sess.valTrials(rwdShift+1)) sess.ind(end)];
    % Epoch into pre/post reward shift for reward shift sessions
    sessFam = epochStruc(sess,1,famInds);
    sessNov = epochStruc(sess,1,novInds);

    % LSI
    [famDat(ct).lckHz, famDat(ct).lapLSI] = get_lickDiscrim(sessFam,[famR1 famR2]*100,zLen);
    [novDat(ct).lckHz, novDat(ct).lapLSI] = get_lickDiscrim(sessNov,[novR1 novR2]*100,zLen);
    % lsiFig = plot_lickDiscrim(sess,[r1pos r2pos]*100,zLen,10);
    uLSIFam(ct) = [mean(famDat(ct).lapLSI(end-29:end))];
    uLSINov(ct) = [mean(novDat(ct).lapLSI(useLaps))];

    % Operant success
    famDat(ct).lapScore = get_trialSuccess(sessFam);
    novDat(ct).lapScore = get_trialSuccess(sessNov);
    % opScoreFig = plot_trialScore(sess);
    ulpSFam(ct) = [mean(famDat(ct).lapScore(end-29:end))];
    ulpSNov(ct) = [mean(novDat(ct).lapScore(useLaps))];

    % % Velocity
    % famDat(ct).vCorXLap = corr(sess.velXpos', mean(sess.velXpos(end-19:end,:))','rows','complete');
    % tmpVCorrFig = figure; hold on; plot(sess.valTrials,bhvr(i).vCorXLap,'k-o',sess.valTrials,smooth(bhvr(i).vCorXLap,5),'r')
    % ylim([-0.5 1]); ylabel('Velocity Corr. to last 20 laps'); xlabel('Lap #'); 
    % set(gca,'FontSize',16,'FontName','Arial')
    % vCorSttEnd(i,:) = [mean(bhvr(i).vCorXLap(lrange(1,:))) mean(bhvr(i).vCorXLap(lrange(2,:)))];

    if saveFlag
        % fsave(lsiFig,[sbase '_lsiXtrial_rzshift_lsi'],1,0,0)
    end
    close all

    ct = ct + 1;
end

cd(sdir)
if saveFlag
    save('AT_bhvr_vars_D5_lsi','famDat','novDat','ulpSFam','ulpSNov','uLSIFam','uLSINov');
end

%% Compare Exp vs Ctl during RZ shift
expMice = logical([1 0 0 1 1]);
vColors = [.35 .35 .35; 1 .25 .25];

tmpLSIFam = plot_barXgroup(uLSIFam(~expMice), uLSIFam(expMice), vColors);
ylim([0 1.1]); ylabel('LSI Familiar')
tmpLSINov = plot_barXgroup(uLSINov(~expMice), uLSINov(expMice), vColors);
ylim([-1.1 0]); ylabel('LSI Novel')

tmpLpsFam = plot_barXgroup(ulpSFam(~expMice), ulpSFam(expMice), vColors);
ylim([0 1.1]); ylabel('Lap Score Familiar')
tmpLpsNov = plot_barXgroup(ulpSNov(~expMice), ulpSNov(expMice), vColors);
ylim([0 1.1]); ylabel('Lap Score Novel')

fsave(tmpLSIFam,[sName 'LSI_fam'],1,0,0)
fsave(tmpLSINov,[sName 'LSI_nov'],1,0,0)
fsave(tmpLpsFam,[sName 'LpS_fam'],1,0,0)
fsave(tmpLpsNov,[sName 'LpS_nov'],1,0,0)

%% Compare across laps and mice

figure; hold on
for i = 1:length(famDat)
    if famDat(i).gExp == 1
        vcol = 'r';
    else
        vcol = 'k';
    end
    plot(smooth(novDat(i).lapLSI),vcol)
end

%% Functions

function [fhandle] = plot_barXgroup(dat1,dat2,vcolors)

nMice = [numel(dat1) numel(dat2)];
uDat = [mean(dat1) mean(dat2)];
stdDat = [std(dat1) std(dat2)];

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.1 0.27])
b = bar(uDat,'FaceColor','flat','HandleVisibility','off');
b.CData = vcolors;
errorbar(1:2,uDat,stdDat ./ sqrt(nMice),'k.','HandleVisibility','off')
plot(1*ones(nMice(1),1),dat1,'ko')
plot(2*ones(nMice(2),1),dat2,'ko')
xlim([0.5 2.5])
xticks(1:2); xticklabels({'Control', 'Exp.'})
set(gca,'FontSize',12,'FontName','Arial')

end

function [fhandle] = plot_bhvrXdays(dat1,dat2,x1,x2,vcolors)

nMice = [numel(dat1) numel(dat2)];
nDays = size(dat1,2);

fhandle = figure; hold on
set(gcf,'units','normalized','Position',[0.3724 0.4056 0.2505 0.3352])
errorbar(1:nDays, mean(dat1), std(dat1) ./ sqrt(nMice(1)),'Linewidth',2,'Color',vcolors(1,:))
errorbar(1:nDays, mean(dat2), std(dat2) ./ sqrt(nMice(2)),'Linewidth',2,'Color',vcolors(2,:))
plot(x1'-0.15,dat1','o','Color',vcolors(1,:))
plot(x2'+0.15,dat2','o','Color',vcolors(2,:))
xlim([0.5 nDays+0.5]); xlabel('Day')
legend('Control','Experimental','Location','best')
set(gca,'FontSize',16,'FontName','Arial')

end