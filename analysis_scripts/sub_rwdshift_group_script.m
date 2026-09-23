%% Prep Sessions/Animals/variables for later blocks
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%

datT = import_xldat("D:\Data\Kelton\analyses\group_analyses","dat_include.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\bigcohort_070926';
cd(groupSDir) 

mInclude = {'KW101','KW097','KW099','KW077','KW073','ZM032','ZM006','ZM035',...
    'KW100','KW082','KW079','ZM012','ZM020','ZM029','KW091','KW087','KW080'}; %By row: Learners, Non-learners
% mInclude = {'KW101','KW097','KW077','KW073','ZM032','ZM006',...
%     'KW100','KW099','HE002','KW082','KW079','KW074','ZM012','ZM029',...
%     'KW094','KW091','KW087','KW080','ZM020'}; % Original split with partials

sessType = 2;
% useInds = datT.include == 1;
for i = 1:height(datT)
    useInds(i) = logical(sum(strcmp(datT.mouse(i),mInclude))) & logical(sum(ismember(sessType,datT.sess_type(i)))); 
end
datT(~useInds,:) = [];  %Clean excluded sessions

saveFlag = 1;
sbase = 'subRwdShift_';
fname = [sbase 'groupDat_phys2'];
bvName = [sbase 'groupDat_bhvr2'];

dbnsz = 0.05;
histoBnsz = 5;
binedges = 0:5:185;
binpos = 0.025:dbnsz:1.825;
wlen = 150;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
fvncols = [.35 .35 .35; 1 .25 .25]; % Gray vs red
% lnlcols = [0 0 1; 0.7656 0.6406 0.5156]; % Blue vs brown
lnlcols = [0.0508 0.4883 0.5273; 0.7852 0.6055 0.2188]; % learn teal vs Nonlearn brown
rgncols = [0.9961 0.7305 0.4336; 0.3672 0.2969 0.3711]; % Sub gold vs CA1 dull purple

clear ps stats

%% Combine Sub RZ Shift Behavior data
combine_bhvrDat(datT,bvName,groupSDir,2);

%% Combine ephys data
combine_rzShiftDat(datT,fname,groupSDir,2,1); % sessType = 2

%% Load previously saved ephys data
cd(groupSDir)
load(fname)

nMice = length(unique(recID(:,1)));
mID = unique(recID(:,1),'stable');
nTotal = size(recID,1);
r1posInd = find(binpos > r1pos,1);
r2posInd = find(binpos > r2pos,1);
sbase = 'subRwdShift_';

% vlFrstID = useCC & vlDat(:,1) <= 0.05;
% vlLastID = useCC & vlDat(:,4) <= 0.05;
% vlBothID = vlFrstID & vlLastID;
% thFrstID = useCC & thDat(:,1) <= 0.05;
% thLastID = useCC & thDat(:,4) <= 0.05;
% thBothID = thFrstID & thLastID;
siFrstID = useCC & lcDat(:,1) <= 0.05;
siLastID = useCC & lcDat(:,5) <= 0.05;
siBothID = siFrstID & siLastID;
siEithID = siFrstID | siLastID;
swrFrstID = useCC & rpDat(:,2) > 1;
swrLastID = useCC & rpDat(:,4) > 1;
swrBothID = swrFrstID & swrLastID;
bstFrstID = useCC & bsDat(:,1) > 0;
bstLastID = useCC & bsDat(:,4) > 0;
bstBothID = bstFrstID & bstLastID;

mLern = [101 99 97 77 73 32 6 35];
mNonl = [100 82 79 12 20 29 91 87 80];
lnInd = [];
nlInd = [];

for i = 1:nMice
    tmpUnits = find(recID(:,1) == mID(i) & rgDat == 2); % rgDat specifies region (1 = CA1; 2 = Sub)
    if ~isempty(find(mLern == mID(i), 1))
        lnInd = [lnInd; tmpUnits];
        mInds(i).grp = 1*ones(size(tmpUnits));
    else
        nlInd = [nlInd; tmpUnits];
        mInds(i).grp = 2*ones(size(tmpUnits));
    end
    mInds(i).logi = false(nTotal,1);
    mInds(i).logi(tmpUnits) = true;
    mInds(i).nSub = sum(recID(:,1) == mID(i) & rgDat == 2 & useCC);
    mInds(i).nCA1 = sum(recID(:,1) == mID(i) & rgDat == 1 & useCC);
end

grp(1).inds = lnInd;
grp(1).mice = mLern;
grp(1).logi = false(nTotal,1);
grp(1).logi(grp(1).inds) = true;
grp(1).sname = '_learn';
grp(2).inds = nlInd;
grp(2).mice = mNonl;
grp(2).logi = false(nTotal,1);
grp(2).logi(grp(2).inds) = true;
grp(2).sname = '_nolearn';

%% Load saved behavior data
cd(groupSDir)
load(bvName)
nMice = length(unique(bhvID(:,1)));
mID = unique(bhvID(:,1),'stable');

mLern = [101 99 97 79 77 73 35 32 6];
mNonl = [100 91 87 82 80 29 20 12];

lnInd = [];
nlInd = [];

for i = 1:size(bhvID,1)
    if ~isempty(find(mLern == bhvID(i,1), 1))
        lnInd = [lnInd; find(bhvID(:,1) == bhvID(i,1))];
    else
        nlInd = [nlInd; find(bhvID(:,1) == bhvID(i,1))];
    end
end

bvgrp(1).bvInd = lnInd;
bvgrp(2).bvInd = nlInd;
bvgrp(1).grpname = 'learn';
bvgrp(2).grpname = 'nolearn';
bvgrp(1).mID = bhvID(bvgrp(1).bvInd);
bvgrp(2).mID = bhvID(bvgrp(2).bvInd);
bvgrp(1).n = length(bvgrp(1).bvInd);
bvgrp(2).n = length(bvgrp(2).bvInd);

%% Behavior comparisons
% Learners: P(Nov Rwd) > 0.55

for i = 1:nMice
    uLapRwd50(i,:) = [mean(bvDat(i).preLapRwd(1:50)) mean(bvDat(i).pstLapRwd(1:50))];
end
% uLDI = [vertcat(bvDat.uPreLckDI), vertcat(bvDat.uPstLckDI)];
uPsv = [vertcat(bvDat.uPreLckPsv), vertcat(bvDat.uPstLckPsv)];
uLapRwd = [vertcat(bvDat.uPreLapRwd), vertcat(bvDat.uPstLapRwd)];
nLaps = [vertcat(bvDat.preNLap), vertcat(bvDat.pstNLap)];

% rwdWtlckSplitF = plot_2d_bhvr(uLDI .* uLapRwd,bvgrp(1).bvInd,[],bvgrp(2).bvInd);
% plot([-1 1],[-0.2 -0.2],'k--')
% plot([0.2 0.2],[-1 1],'k--')
% xlabel('uFam P(Rwd)*LSI'); xlim([-1 1])
% ylabel('uNov P(Rwd)*LSI'); ylim([-1 1])
% legend('Learner','Non-learner')
%
% ldiSplitF = plot_2d_bhvr(uLDI,bvgrp(1).bvInd,[],bvgrp(2).bvInd);
% plot([-1 1],[-0.4 -0.4],'k--')
% xlabel('mean Fam LSI'); xlim([-1 1])
% ylabel('mean Nov LSI'); ylim([-1 1])
% legend('Learner','Non-learner')
% 
% combiSplitF = figure;
% plot3(uLDI(bvgrp(1).bvInd,2),nLaps(bvgrp(1).bvInd,2),uLapRwd(bvgrp(1).bvInd,2),'bo',...
%     uLDI(bvgrp(2).bvInd,2),nLaps(bvgrp(2).bvInd,2),uLapRwd(bvgrp(2).bvInd,2),'ro');
% xlabel('LSI Nov');
% ylabel('N Laps Nov');
% zlabel('P(Nov Lap Rewarded');

psvSplitF = plot_2d_bhvr(uPsv,bvgrp(1).bvInd,bvgrp(2).bvInd,lnlcols);
% plot([-1 1],[-0.15 -0.15],'k--')
xlabel('LDI F'); xlim([-1 1])
ylabel('LDI N'); ylim([-1 1])

rwdSplitF = plot_2d_bhvr(uLapRwd,bvgrp(1).bvInd,bvgrp(2).bvInd,lnlcols);
plot([0 1],[0.55 0.55],'k--')
xlabel('P(F Lap Rewarded)'); xlim([0 1])
ylabel('P(N Lap Rewarded)'); ylim([0 1])
legend({'Strong adapter','Weak adapter'}, 'Location','sw')

lapSplitF = plot_2d_bhvr(nLaps,bvgrp(1).bvInd,bvgrp(2).bvInd,lnlcols);
xlabel('# Laps F'); xlim([0 140])
ylabel('# Laps N'); ylim([0 140])

if saveFlag
    % fsave(rwdWtlckSplitF,[sbase 'bhv_wtlckPrePst'],1,0);
    % fsave(ldiSplitF,[sbase 'bhv_lckPrePst'],1,0);
    % fsave(combiSplitF,[sbase 'bhv_combiPst'],1,0);
    fsave(psvSplitF,[sbase 'bhv_psvPrePst'],1,1);
    fsave(rwdSplitF,[sbase 'bhv_rwdPrePst'],1,1);
    fsave(lapSplitF,[sbase 'bhv_lapPrePst'],1,1);
end

%% Statistical Quantification - Learn Vs Non Learn

[~,ps.ttest_psv_fam,~,stats.ttest_psv_fam] = ttest2(uPsv(bvgrp(1).bvInd,1),uPsv(bvgrp(2).bvInd,1));
[~,ps.ttest_psv_nov,~,stats.ttest_psv_nov] = ttest2(uPsv(bvgrp(1).bvInd,2),uPsv(bvgrp(2).bvInd,2));
[~,ps.ttest_lap_fam,~,stats.ttest_lap_fam] = ttest2(nLaps(bvgrp(1).bvInd,1),nLaps(bvgrp(2).bvInd,1));
[~,ps.ttest_lap_nov,~,stats.ttest_lap_nov] = ttest2(nLaps(bvgrp(1).bvInd,2),nLaps(bvgrp(2).bvInd,2));
[~,ps.ttest_rwd_fam,~,stats.ttest_rwd_fam] = ttest2(uLapRwd(bvgrp(1).bvInd,1),uLapRwd(bvgrp(2).bvInd,1));
[~,ps.ttest_rwd_nov,~,stats.ttest_rwd_nov] = ttest2(uLapRwd(bvgrp(1).bvInd,2),uLapRwd(bvgrp(2).bvInd,2));

psvScatF = plot_grpBhv(uPsv,bvgrp,lnlcols);
xticklabels({'F','N'}); ylim([-1 1]); 
text2bar(psvScatF,'',ps.ttest_psv_fam,0.85, 0.5); text2bar(psvScatF,'LDI',ps.ttest_psv_nov,0.25);

rwdScatF = plot_grpBhv(uLapRwd,bvgrp,lnlcols);
xticklabels({'F','N'}); ylim([0 1.1]); 
text2bar(rwdScatF,'',ps.ttest_rwd_fam,0.85,0.5); text2bar(rwdScatF,'P(Rwd)',ps.ttest_rwd_nov,0.25);

lapScatF = plot_grpBhv(nLaps,bvgrp,lnlcols);
xticklabels({'F','N'}); ylim([0 140]); 
text2bar(lapScatF,'',ps.ttest_lap_fam,0.85); text2bar(lapScatF,'# Laps',ps.ttest_lap_nov,0.25);

% lapBarF = plotMiniBar(nLaps(bvgrp(1).bvInd,1),nLaps(bvgrp(2).bvInd,1),vColors);
% xticklabels({'Learn','Non-Lrn'}); text2bar(lapBarF,'mean Fam Laps',ps.ttest_lap_fam); ylim([0 170])
% lapBarN = plotMiniBar(nLaps(bvgrp(1).bvInd,2),nLaps(bvgrp(2).bvInd,2),vColors);
% xticklabels({'Learn','Non-Lrn'}); text2bar(lapBarN,'mean Nov Laps',ps.ttest_lap_nov); ylim([0 170])
% rwdBarF = plotMiniBar(uLapRwd(bvgrp(1).bvInd,1),uLapRwd(bvgrp(2).bvInd,1),vColors);
% xticklabels({'Learn','Non-Lrn'}); text2bar(rwdBarF,'mean Fam P(rwd)',ps.ttest_rwd_fam); ylim([0 1])
% rwdBarN = plotMiniBar(uLapRwd(bvgrp(1).bvInd,2),uLapRwd(bvgrp(2).bvInd,2),vColors);
% xticklabels({'Learn','Non-Lrn'}); text2bar(rwdBarN,'mean Nov P(rwd)',ps.ttest_rwd_nov); ylim([0 1])
% psvBarF = plotMiniBar(uPsv(bvgrp(1).bvInd,1),uPsv(bvgrp(2).bvInd,1),vColors);
% xticklabels({'Learn','Non-Lrn'}); text2bar(psvBarF,'mean Fam Perseveration',ps.ttest_psv_fam); ylim([-1 1])
% psvBarN = plotMiniBar(uPsv(bvgrp(1).bvInd,2),uPsv(bvgrp(2).bvInd,2),vColors);
% xticklabels({'Learn','Non-Lrn'}); text2bar(psvBarN,'mean Nov Perseveration',ps.ttest_psv_nov); ylim([-1 1])

if saveFlag
    fsave(psvScatF,[sbase 'bhv_psv_errbar'],1,1)
    fsave(rwdScatF,[sbase 'bhv_rwd_errbar'],1,1)
    fsave(lapScatF,[sbase 'bhv_lap_errbar'],1,1)
end

%% Behavior comparisons first L1-10 vs L41-50 in F vs N
xsF = 1:50;
xsN = 51:100;

for i = 1:length(lnInd)
    % trialLDIF(i,1:nLaps(lnInd(i),1)) = bvDat(lnInd(i)).preLckDI;
    % trialLDIN(i,1:nLaps(lnInd(i),2)) = bvDat(lnInd(i)).pstLckDI;
    trialPsvF_ln(i,1:nLaps(lnInd(i),1)) = bvDat(lnInd(i)).preLckPsv;
    trialPsvN_ln(i,1:nLaps(lnInd(i),2)) = bvDat(lnInd(i)).pstLckPsv;
end
for i = 1:length(nlInd)
    trialPsvF_nl(i,1:nLaps(nlInd(i),1)) = bvDat(nlInd(i)).preLckPsv;
    trialPsvN_nl(i,1:nLaps(nlInd(i),2)) = bvDat(nlInd(i)).pstLckPsv;
end

mdlFPsv_ln = get_linfit(xsF,mean(trialPsvF_ln(:,xsF)));
mdlNPsv_ln = get_linfit(xsF,mean(trialPsvN_ln(:,xsF)));
mdlFPsv_nl = get_linfit(xsF,mean(trialPsvF_nl(:,xsF)));
mdlNPsv_nl = get_linfit(xsF,mean(trialPsvN_nl(:,xsF)));

[ciup_F_ln, cidn_F_ln] = get_CI(trialPsvF_ln(:,xsF));
[ciup_N_ln, cidn_N_ln] = get_CI(trialPsvN_ln(:,xsF));
[ciup_F_nl, cidn_F_nl] = get_CI(trialPsvF_nl(:,xsF));
[ciup_N_nl, cidn_N_nl] = get_CI(trialPsvN_nl(:,xsF));

% mdlF = get_linfit(1:50,mean(trialLDIF(:,1:50)));
% mdlN = get_linfit(1:50,mean(trialLDIN(:,1:50)));
% 
% corLapLDIFNF = figure; hold on
% plot(1:50, mean(trialLDIF(:,1:50)),'k', 51:100, mean(trialLDIN(:,1:50)),'r')
% plot(1:50,mdlF.ypred, 'k', 'LineWidth',2)
% plot(51:100,mdlN.ypred,'k','LineWidth',2)
% xlabel("Lap #")
% text2bar(corLapLDIFNF,"Mean LDI",mdlF.p,0.8,0.4);
% text2bar(corLapLDIFNF,"Mean LDI",mdlN.p,0.3,0.5);
% set(gca,'FontSize',12,'FontName','Arial')

corLapPsvLnF = figure; hold on
set(gcf,'Units','normalized','Position',[1.2 0.4 0.6315 0.1703])
plot_CIs(xsF, ciup_F_ln, cidn_F_ln, lnlcols(1,:)/2);
plot_CIs(xsN, ciup_N_ln, cidn_N_ln, lnlcols(1,:));
plot(xsF, mean(trialPsvF_ln(:,1:50)), 'color', lnlcols(1,:)/2)
plot(xsN, mean(trialPsvN_ln(:,1:50)), 'color', lnlcols(1,:))
plot(xsF,mdlFPsv_ln.ypred, 'color', lnlcols(1,:)/2, 'LineWidth',2)
plot(xsN,mdlNPsv_ln.ypred, 'color', lnlcols(1,:),   'LineWidth',2)
xlabel("Lap #"); ylim([-1.1 1.1])
text2bar(corLapPsvLnF,"",    mdlFPsv_ln.p, 0.8, 0.3, lnlcols(1,:)/2);
text2bar(corLapPsvLnF,"LDI", mdlNPsv_ln.p, 0.3, 0.8, lnlcols(1,:));
set(gca,'FontSize',16,'FontName','Arial')

corLapPsvNlF = figure; hold on
set(gcf,'Units','normalized','Position',[1.2 0.4 0.6315 0.1703])
plot_CIs(xsF, ciup_F_nl, cidn_F_nl, lnlcols(2,:)/2);
plot_CIs(xsN, ciup_N_nl, cidn_N_nl, lnlcols(2,:));
plot(xsF, mean(trialPsvF_nl(:,1:50)), 'color', lnlcols(2,:)/2)
plot(xsN, mean(trialPsvN_nl(:,1:50)), 'color', lnlcols(2,:))
plot(xsF,mdlFPsv_nl.ypred, 'color', lnlcols(2,:)/2, 'LineWidth',2)
plot(xsN,mdlNPsv_nl.ypred, 'color', lnlcols(2,:),   'LineWidth',2)
xlabel("Lap #"); ylim([-1.1 1.1])
text2bar(corLapPsvNlF,"",    mdlFPsv_nl.p, 0.8, 0.5, lnlcols(2,:)/2);
text2bar(corLapPsvNlF,"LDI", mdlNPsv_nl.p, 0.3, 0.2, lnlcols(2,:));
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    % fsave(corLapLDIFNF,[sbase 'bhv_lapXldi_corr'],1,0)
    fsave(corLapPsvLnF,[sbase 'bhv_lapXpsv_corr_learn'],1,1)
    fsave(corLapPsvNlF,[sbase 'bhv_lapXpsv_corr_nolrn'],1,1)
end

%% Averaged lick and velocity profiles
lBins = 0:0.03:1.85;
mBin = round(length(lBins)/2);

for i = 1:nMice
    preRZBin = find(lBins > bvDat(i).rzPos(1),1);
    pstRZBin = find(lBins > bvDat(i).rzPos(2),1);
    
    uPreVel(i,:) = mean(bvDat(i).preSpVelMap,'omitnan');
    uPreVel(i,:) = circshift(uPreVel(i,:), mBin - preRZBin);
    uPstVel(i,:) = mean(bvDat(i).pstSpVelMap,'omitnan');
    uPstVel(i,:) = circshift(uPstVel(i,:), mBin - pstRZBin);
    
    uPreLck(i,:) = mean(bvDat(i).preLMap,'omitnan');
    uPreLck(i,:) = circshift(uPreLck(i,:), mBin - preRZBin);
    uPstLck(i,:) = mean(bvDat(i).pstLMap,'omitnan');
    uPstLck(i,:) = circshift(uPstLck(i,:), mBin - pstRZBin);
end

% Average group types pre and post separately
vColorsLN = [0 0 1; 0.7656 0.6406 0.5156];
vColorsFN = [0 0 0; 1 0 0];

groupLckPreF = plot_3bhvrTraceCI(uPreLck,bvgrp,vColorsLN);
ylim([0 15]); ylabel('Licks (Hz)'); legend('Learner','Non-learner')
groupLckPstF = plot_3bhvrTraceCI(uPstLck,bvgrp,vColorsLN);
ylim([0 15]); ylabel('Licks (Hz)')

groupVelPreF = plot_3bhvrTraceCI(uPreVel,bvgrp,vColorsLN);
xlim([-100 100]); ylim([0 45]); ylabel('Velocity (cm/s)')
groupVelPstF = plot_3bhvrTraceCI(uPstVel,bvgrp,vColorsLN);
xlim([-100 100]); ylim([0 45]); ylabel('Velocity (cm/s)')

% Plot pre/post by group
uLckLernF = plot_bhvrTraceCI(uPreLck(lnInd,:),uPstLck(lnInd,:),vColorsFN);
ylim([0 15]); ylabel('Licks (Hz)'); legend('Familiar','Novel')
uLckNonLF = plot_bhvrTraceCI(uPreLck(nlInd,:),uPstLck(nlInd,:),vColorsFN);
ylim([0 15]); ylabel('Licks (Hz)')

uVelLernF = plot_bhvrTraceCI(uPreVel(lnInd,:),uPstVel(lnInd,:),vColorsFN);
ylim([0 45]); ylabel('Velocity (cm/s)')
uVelNonLF = plot_bhvrTraceCI(uPreVel(nlInd,:),uPstVel(nlInd,:),vColorsFN);
ylim([0 45]); ylabel('Velocity (cm/s)')

if saveFlag
    fsave(groupLckPreF,[sbase 'bhv_meanLckXgrp_pre'],1,1)
    fsave(groupLckPstF,[sbase 'bhv_meanLckXgrp_pst'],1,1)
    fsave(groupVelPreF,[sbase 'bhv_meanVelXgrp_pre'],1,1)
    fsave(groupVelPstF,[sbase 'bhv_meanVelXgrp_pst'],1,1)
    fsave(uLckLernF,[sbase 'bhv_meanLck_prepst_learn'],1,1)
    fsave(uLckNonLF,[sbase 'bhv_meanLck_prepst_nonlearn'],1,1)
    fsave(uVelLernF,[sbase 'bhv_meanVel_prepst_learn'],1,1)
    fsave(uVelNonLF,[sbase 'bhv_meanVel_prepst_nonlearn'],1,1)
end

%% Statistical comparison behavior pre vs post shift

for i = 1:2
    [~,ps(i).ttest_psv_prepst,~,stats(i).ttest_psv_prepst] = ttest(uPsv(bvgrp(i).bvInd,1), uPsv(bvgrp(i).bvInd,2));
    [~,ps(i).ttest_rwd_prepst,~,stats(i).ttest_rwd_prepst] = ttest(uLapRwd(bvgrp(i).bvInd,1), uLapRwd(bvgrp(i).bvInd,2));
    [~,ps(i).ttest_lap_prepst,~,stats(i).ttest_lap_prepst] = ttest(nLaps(bvgrp(i).bvInd,1), nLaps(bvgrp(i).bvInd,2));

    psvPrePstF = plot_barXmouse(uPsv(bvgrp(i).bvInd,:),fvncols);
    text2bar(psvPrePstF,'LDI',ps(i).ttest_psv_prepst); legend('off')

    rwdPrePstF = plot_barXmouse(uLapRwd(bvgrp(i).bvInd,:),fvncols);
    text2bar(rwdPrePstF,'P(Operant Reward)',ps(i).ttest_rwd_prepst); legend('off')

    lapPrePstF = plot_barXmouse(nLaps(bvgrp(i).bvInd,:),fvncols);
    text2bar(lapPrePstF,'# Laps',ps(i).ttest_lap_prepst);

    if saveFlag
        fsave(psvPrePstF,[sbase 'bhv_meanPsv_prepst_' bvgrp(i).grpname],1,1)
        fsave(rwdPrePstF,[sbase 'bhv_meanRwd_prepst_' bvgrp(i).grpname],1,1)
        fsave(lapPrePstF,[sbase 'bhv_meanLap_prepst_' bvgrp(i).grpname],1,1)
        close all
    end
end

%% Start Ephys comparisons
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================

%% Compare unit counts Learn vs Non learn

nCCmat = [vertcat(mInds.nSub), vertcat(mInds.nCA1)];
nCCmat(nCCmat == 0) = nan;

[~,ps(i).nCC_sub_lnnl,~,stats(i).nCC_sub_lnnl] = ttest2(nCCmat(bvgrp(1).bvInd,1), nCCmat(bvgrp(2).bvInd,1));
[~,ps(i).nCC_ca1_lnnl,~,stats(i).nCC_ca1_lnnl] = ttest2(nCCmat(bvgrp(1).bvInd,2), nCCmat(bvgrp(2).bvInd,2));

nSubCompF = plotBar2(nCCmat(bvgrp(1).bvInd,1), nCCmat(bvgrp(2).bvInd,1),lnlcols);
text2bar(nSubCompF,'# of Sub units',ps(i).nCC_sub_lnnl); ylim([0 150]); xticklabels({'Strong', 'Weak'})
nCA1CompF = plotBar2(nCCmat(bvgrp(1).bvInd,2), nCCmat(bvgrp(2).bvInd,2),lnlcols);
text2bar(nCA1CompF,'# of CA1 units',ps(i).nCC_ca1_lnnl); ylim([0 150]); xticklabels({'Strong', 'Weak'})

if saveFlag
    fsave(nSubCompF, [sbase 'nCC_sub_lnnl']);
    fsave(nCA1CompF, [sbase 'nCC_ca1_lnnl']);
end

%% Compare FR, SI, and Peak rate within each behavior group
% lcDat: 1&5 = sig.; 2&6 = SI; 3&7 = pkRate; 4&8 = pkLoc
% frDat: 1&3 = standing; 2&4 = running

for i = 1:2
    % Test FR for all in-layer pyramidal units
    [~,ps(i).fr_PPStnd_all,~,stats(i).fr_PPStnd_all] = ttest(frDat(grp(i).logi & useCC,1),frDat(grp(i).logi & useCC,3));
    [~,ps(i).fr_PPRunn_all,~,stats(i).fr_PPRunn_all] = ttest(frDat(grp(i).logi & useCC,2),frDat(grp(i).logi & useCC,4));

    fnames(i).frStndRunFig = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
    fnames(i).frStndRunFig = fixRatio(fnames(i).frStndRunFig);
    violinplot(frDat(grp(i).logi & useCC,[1 3 2 4]), ones(sum(grp(i).logi & useCC),1), 'ViolinColor',[fvncols(1,:) / 2; fvncols(2,:) / 2; fvncols(1,:); fvncols(2,:)],'ShowData',false);
    xlim([0.5 4.5]); xticklabels({'F Stand','N Stand','F Run','N Run'})
    text2bar(fnames(i).frStndRunFig,'Firing Rate (Hz)',ps(i).fr_PPStnd_all,0.85);
    text2bar(fnames(i).frStndRunFig,'Firing Rate (Hz)',ps(i).fr_PPRunn_all,0.35);
    set(gca,'FontSize',12,'FontName','Arial')

    % Quantify proportions of spatially modulated units
    [fnames(i).siPie, fnames(i).siprcts] = prepostPie(siFrstID(grp(i).inds),siLastID(grp(i).inds),useCC(grp(i).inds));
    title("Sig. Spatial Info. units");
    siBothCounts = [groupcounts(recID(siFrstID & grp(i).logi,1),sort(grp(i).mice),'IncludeEmptyGroups',true),...
                    groupcounts(recID(siLastID & grp(i).logi,1),sort(grp(i).mice),'IncludeEmptyGroups',true)];
    siBothRatio = siBothCounts ./ [groupcounts(recID(useCC & grp(i).logi,1), sort(grp(i).mice),'IncludeEmptyGroups',true) groupcounts(recID(useCC & grp(i).logi,1),sort(grp(i).mice),'IncludeEmptyGroups',true)];
    [~,ps(i).lc_SICt_both,~,stats(i).lc_SICt_both] = ttest(siBothCounts(:,1),siBothCounts(:,2));
    [~,ps(i).lc_SICt_both,~,stats(i).lc_SICt_both] = ttest(siBothCounts(:,1),siBothCounts(:,2));
    fnames(i).lcBothCtFig = plot_barXmouse(siBothRatio);
    text2bar(fnames(i).lcBothCtFig,'P(Sig. Spatial Info.)',ps(i).lc_SICt_both);

    % Test SI for units modulated in both phases
    [~,ps(i).lc_FNSI_both,~,stats(i).lc_FNSI_both] = ttest(lcDat(siBothID & grp(i).logi,2), lcDat(siBothID & grp(i).logi,6));
    [~,ps(i).lc_FNPk_both,~,stats(i).lc_FNPk_both] = ttest(lcDat(siBothID & grp(i).logi,3), lcDat(siBothID & grp(i).logi,7));

    fnames(i).lcBothSIFig = plotBar2(lcDat(siBothID  & grp(i).logi,2),lcDat(siBothID  & grp(i).logi,6)); ylim([0 4]);
    text2bar(fnames(i).lcBothSIFig,'Spatial Information (Bits/spike)',ps(i).lc_FNSI_both);
    fnames(i).lcBothPkFig = plotBar2(lcDat(siBothID  & grp(i).logi,3),lcDat(siBothID  & grp(i).logi,7));
    text2bar(fnames(i).lcBothPkFig,'Peak Field FR (Hz)',ps(i).lc_FNPk_both);

    if saveFlag
        fsave(fnames(i).siPie,[sbase 'si_Mod_pie' grp(i).sname])
        fsave(fnames(i).lcBothSIFig,[sbase 'si_both_bar' grp(i).sname])
        fsave(fnames(i).lcBothPkFig,[sbase 'pk_both_bar' grp(i).sname])
        fsave(fnames(i).frStndRunFig,[sbase 'fr_standrun_violin' grp(i).sname])
        % fsave(fnames(i).lcEithSIFig,[sbase 'si_eith_bar' grp(i).sname])
        % fsave(fnames(i).lcEithPkFig,[sbase 'pk_eith_bar' grp(i).sname])
        fsave(fnames(i).lcBothCtFig,[sbase 'si_Mod_both_bar' grp(i).sname])
        close all
    end
end

%% Test FR by region and adapter group, standing/running with LME

frMat = nan(nMice,8);  % Sub F stand, F Run, N Stand, N Run, CA1 F stand, F Run, N Stand, N Run
for i = 1:nMice
    useSub = recID(:,1) == mID(i) & rgDat == 2 & useCC;
    useCA1 = recID(:,1) == mID(i) & rgDat == 1 & useCC;
    frMat(i,:) = [mean(frDat(useSub,:)) mean(frDat(useCA1,:))];
end

frVarNames = {'fr','epoch','grp','mouse'};
[fr_subFNstd_lme_table] = get_lmetable([frMat(:,1); frMat(:,3)],bvgrp,mID,frVarNames);
fr_subFNstd_lme = fitlme(fr_subFNstd_lme_table,'fr ~ epoch * grp + (1|mouse)');
fr_subFNstd_lmeF = plot_2wayLME(frMat(:,1)',frMat(:,3)',bvgrp,lnlcols); ylim([0 20])
xticklabels({'F','N','F','N'})

[fr_subFNrun_lme_table] = get_lmetable([frMat(:,2); frMat(:,4)],bvgrp,mID,frVarNames);
fr_subFNrun_lme = fitlme(fr_subFNrun_lme_table,'fr ~ epoch * grp + (1|mouse)');
fr_subFNrun_lmeF = plot_2wayLME(frMat(:,2)',frMat(:,4)',bvgrp,lnlcols); ylim([0 20])
xticklabels({'F','N','F','N'})

[fr_ca1FNstd_lme_table] = get_lmetable([frMat(:,5); frMat(:,7)],bvgrp,mID,frVarNames);
fr_ca1FNstd_lme = fitlme(fr_ca1FNstd_lme_table,'fr ~ epoch * grp + (1|mouse)');
fr_ca1FNstd_lmeF = plot_2wayLME(frMat(:,5)',frMat(:,7)',bvgrp,lnlcols); ylim([0 20])
xticklabels({'F','N','F','N'})

[fr_ca1FNrun_lme_table] = get_lmetable([frMat(:,6); frMat(:,8)],bvgrp,mID,frVarNames);
fr_ca1FNrun_lme = fitlme(fr_ca1FNrun_lme_table,'fr ~ epoch * grp + (1|mouse)');
fr_ca1FNrun_lmeF = plot_2wayLME(frMat(:,6)',frMat(:,8)',bvgrp,lnlcols); ylim([0 20])
xticklabels({'F','N','F','N'})

% Test avg FR subiculum vs CA1 over all mice
[~,ps.fr_rgn,~,stats.fr_rgn] = ttest2(mean(frMat(:,5:8),2),mean(frMat(:,1:4),2));
fr_rgn_F = plotBar2(mean(frMat(:,5:8),2),mean(frMat(:,1:4),2),rgncols); ylim([0 20])
text2bar(fr_rgn_F,'Global FR (Hz)',ps.fr_rgn); xticklabels({'CA1','Sub'})

if saveFlag
    fsave(fr_subFNstd_lmeF,[sbase 'fr_sub_lme_stand'])
    fsave(fr_subFNrun_lmeF,[sbase 'fr_sub_lme_run'])
    fsave(fr_ca1FNstd_lmeF,[sbase 'fr_ca1_lme_stand'])
    fsave(fr_ca1FNrun_lmeF,[sbase 'fr_ca1_lme_run'])
    fsave(fr_rgn_F,[sbase 'fr_rgn_bar'])
end

%% Compare FR, SI and Peak rate across each behavior group

lnlcols = [0 0 1; 0.7656 0.6406 0.5156];

[~,ps(1).lc_LnNlSI_F_both,~,stats(1).lc_LnNlSI_F_both] = ttest2(lcDat(siBothID & grp(1).logi,2), lcDat(siBothID & grp(2).logi,2));
[~,ps(1).lc_LnNlSI_N_both,~,stats(1).lc_LnNlSI_N_both] = ttest2(lcDat(siBothID & grp(1).logi,6), lcDat(siBothID & grp(2).logi,6));
[~,ps(1).lc_LnNlSI_F_frst,~,stats(1).lc_LnNlSI_F_frst] = ttest2(lcDat(siFrstID & grp(1).logi,2), lcDat(siFrstID & grp(2).logi,2));
[~,ps(1).lc_LnNlSI_N_last,~,stats(1).lc_LnNlSI_N_last] = ttest2(lcDat(siLastID & grp(1).logi,6), lcDat(siLastID & grp(2).logi,6));
[~,ps(1).fr_LnNl_stnd_F,~,stats(1).fr_LnNl_stnd_F] = ttest2(frDat(useCC & grp(1).logi,1), frDat(useCC & grp(2).logi,1));
[~,ps(1).fr_LnNl_runn_F,~,stats(1).fr_LnNl_runn_F] = ttest2(frDat(useCC & grp(1).logi,3), frDat(useCC & grp(2).logi,3));
[~,ps(1).fr_LnNl_stnd_N,~,stats(1).fr_LnNl_stnd_N] = ttest2(frDat(useCC & grp(1).logi,2), frDat(useCC & grp(2).logi,2));
[~,ps(1).fr_LnNl_runn_N,~,stats(1).fr_LnNl_runn_N] = ttest2(frDat(useCC & grp(1).logi,4), frDat(useCC & grp(2).logi,4));

lcBothSI_LnNl_F = plotBar2(lcDat(siBothID & grp(1).logi,2), lcDat(siBothID & grp(2).logi,2),lnlcols); ylim([0 4])
xticklabels({'Strong','Weak'}); text2bar(lcBothSI_LnNl_F,'Spatial Information (Bits/spike)',ps(1).lc_LnNlSI_F_both);
lcBothSI_LnNl_N = plotBar2(lcDat(siBothID & grp(1).logi,6), lcDat(siBothID & grp(2).logi,6),lnlcols); ylim([0 4])
xticklabels({'Strong','Weak'}); text2bar(lcBothSI_LnNl_N,'Spatial Information (Bits/spike)',ps(1).lc_LnNlSI_N_both);
lcFrstSI_LnNl_F = plotBar2(lcDat(siFrstID & grp(1).logi,2), lcDat(siFrstID & grp(2).logi,2),lnlcols); ylim([0 4])
xticklabels({'Strong','Weak'}); text2bar(lcFrstSI_LnNl_F,'Spatial Information (Bits/spike)',ps(1).lc_LnNlSI_F_frst);
lcLastSI_LnNl_N = plotBar2(lcDat(siLastID & grp(1).logi,6), lcDat(siLastID & grp(2).logi,6),lnlcols); ylim([0 4])
xticklabels({'Strong','Weak'}); text2bar(lcLastSI_LnNl_N,'Spatial Information (Bits/spike)',ps(1).lc_LnNlSI_N_last);

if saveFlag
    fsave(lcBothSI_LnNl_F,[sbase 'si_bothF_bar_LnNl'])
    fsave(lcBothSI_LnNl_N,[sbase 'si_bothN_bar_LnNl'])
    fsave(lcFrstSI_LnNl_F,[sbase 'si_frstF_bar_LnNl'])
    fsave(lcLastSI_LnNl_N,[sbase 'si_lastN_bar_LnNl'])
end

%% Burst Metrics pre/post and Sub vs CA1
% bsDat: 1:3 Fam burstIndex, burstISI, burstLen, 4:6 Novel

bsMat = nan(nMice,12);  % Sub F, N, CA1 F, N
for i = 1:nMice
    useSub = recID(:,1) == mID(i) & rgDat == 2 & useCC; % Only Sub Pyrs
    useCA1 = recID(:,1) == mID(i) & rgDat == 1 & useCC;
    bsMat(i,:) = [mean(bsDat(useSub,:),'omitnan') mean(bsDat(useCA1,:),'omitnan')];
end

% Test avg burst rate subiculum vs CA1 over all mice
[~,ps.bs_rgn,~,stats.bs_rgn] = ttest2(mean(bsMat(:,[7 10]),2), mean(bsMat(:,[1 4]),2));
bst_rgn_F = plotBar2(mean(bsMat(:,[7 10]),2), mean(bsMat(:,[1 4]),2), rgncols);
text2bar(bst_rgn_F,'Burst Index',ps.bs_rgn); xticklabels({'CA1','Sub'})

% Sub LMEs for burst Index, ISI, and burstlen between learn / non learners
bsVarNames = {'burstMetric','epoch','grp','mouse'};
[bst_subFNBI_lme_table] = get_lmetable([bsMat(:,1); bsMat(:,4)],bvgrp,mID,bsVarNames);
bst_subFNBI_lme = fitlme(bst_subFNBI_lme_table,'burstMetric ~ epoch * grp + (1|mouse)');
bst_subFNBI_lmeF = plot_2wayLME(bsMat(:,1)', bsMat(:,4)',bvgrp,lnlcols); ylim([0 1])
xticklabels({'F','N','F','N'}); ylabel('Burst Index')

[bst_subFNisi_lme_table] = get_lmetable([bsMat(:,2); bsMat(:,5)],bvgrp,mID,bsVarNames);
bst_subFNisi_lme = fitlme(bst_subFNisi_lme_table,'burstMetric ~ epoch * grp + (1|mouse)');
bst_subFNisi_lmeF = plot_2wayLME(bsMat(:,2)', bsMat(:,5)',bvgrp,lnlcols); ylim([0 7])
xticklabels({'F','N','F','N'}); ylabel('Burst ISI (ms)')

[bst_subFNLen_lme_table] = get_lmetable([bsMat(:,3); bsMat(:,6)],bvgrp,mID,bsVarNames);
bst_subFNLen_lme = fitlme(bst_subFNLen_lme_table,'burstMetric ~ epoch * grp + (1|mouse)');
bst_subFNLen_lmeF = plot_2wayLME(bsMat(:,3)', bsMat(:,6)',bvgrp,lnlcols); ylim([0 7])
xticklabels({'F','N','F','N'}); ylabel('Burst Len (spikes)')

if saveFlag
    fsave(bst_rgn_F,[sbase 'bst_rgn_bar'])
    fsave(bst_subFNBI_lmeF,[sbase 'bst_sub_lme_burstInd'])
    fsave(bst_subFNisi_lmeF,[sbase 'bst_sub_lme_burstISI'])
    fsave(bst_subFNLen_lmeF,[sbase 'bst_sub_lme_burstLen'])
end

%% Waterfall by group RZ at 0

unitgrp = siBothID; % Inclusion mask based on SI

if unitgrp == siBothID
    grptag = 'both';
elseif unitgrp == siEithID
    grptag = 'eith';
end

mBin = round(length(binpos)/2);
subGrp = [];
binedges = -0.90:0.05:0.95;
nvals = 3;  % 3xbinsize for peak distribution
binedges2 = binedges(1:nvals :end); % For peak distribution plots
nBins = length(binedges) - 1;
xcoords = (binedges2(1:end-1) + 0.5*diff(binedges2(1:2)))*100;
rShiftUnit = zeros(length(useCC),1);

clear mousePkPreSum mousePkPstSum
for i = 1:nMice
    preRZBin = find(binpos > bvDat(i).rzPos(1),1);
    pstRZBin = find(binpos > bvDat(i).rzPos(2),1);
    tmpUnits = unitgrp & mInds(i).logi; % Only pyrs in region & SI unit group
    try
        subGrp  = [subGrp; ones(sum(tmpUnits),1)*mInds(i).grp(1)];
    catch
    end
    rShiftUnit(tmpUnits) = mBin - pstRZBin;
    alnMouseMap(i).pre = circshift(lcMap(tmpUnits,1:length(binpos)),     mBin - preRZBin, 2);
    alnMouseMap(i).pst = circshift(lcMap(tmpUnits,length(binpos)+1:end), mBin - pstRZBin, 2);
    [~,~,mousePkPre(i,:)] = plot_unitPkHisto(alnMouseMap(i).pre,binedges2*100,1,0);
    [~,~,mousePkPst(i,:)] = plot_unitPkHisto(alnMouseMap(i).pst,binedges2*100,1,0);
    mousePkPreSum(i,:) = arrayfun(@(x) sum(mousePkPre(i,x:x+nvals-1)),1:nvals:size(mousePkPre,2)-nvals+1); % Combine bins for cleaner counts
    mousePkPstSum(i,:) = arrayfun(@(x) sum(mousePkPst(i,x:x+nvals-1)),1:nvals:size(mousePkPst,2)-nvals+1);
end
alnSpMapPre = vertcat(alnMouseMap.pre);
alnSpMapPst = vertcat(alnMouseMap.pst);

for i = 1:2
    alnSpMap(i).pre = alnSpMapPre(subGrp == i,:);
    alnSpMap(i).pst = alnSpMapPst(subGrp == i,:);

    [fnames(i).spBothPreSortPreFig,tmpMap,mSort(i).sortPre] = plot_unitWaterfall(alnSpMap(i).pre,binedges,0,1,0);
    % xticks([1, round(nBins)/2, nBins]); xticklabels([-90, 0, 90])
    plot([mBin mBin],[0 sum(unitgrp)],'k--','LineWidth',2); title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')

    % [fnames(i).spBothPreSortPreHisto,pkMapPre] = plot_unitPkHisto(tmpMap,binedges*100,1,0);
    % xticks([binedges(1),binedges(round(nBins/2)),binedges(end)]*100); xticklabels([-90, 0, 90])
    % plot([0 0],[0 0.11],'k--','LineWidth',2); xlabel('Track Position (cm)'); xlim([-90 95]); ylim([0 0.1])
    % [ps(i).lc_bothPrePkUniformity, stats.lc_bothPrePkUniformity] = pkChi2(pkMapPre,binedges);
    % text2bar(fnames(i).spBothPreSortPreHisto,'',ps(i).lc_bothPrePkUniformity,0.9)

    [fnames(i).spBothPstSortPstFig,tmpMap] = plot_unitWaterfall(alnSpMap(i).pst,binedges,0,1,0);
    plot([mBin mBin],[0 sum(unitgrp)],'r--','LineWidth',2); title('Novel RZ, sort Novel'); xlabel('Track Position (cm)')

    % [fnames(i).spBothPstSortPstHisto,pkMapPst] = plot_unitPkHisto(tmpMap,binedges*100,1,0);
    % plot([0 0],[0 0.11],'r--','LineWidth',2); xlabel('Track Position (cm)'); xlim([-90 95]); ylim([0 0.1])
    % [ps(i).lc_bothPstPkUniformity, stats.lc_bothPstPkUniformity] = pkChi2(pkMapPst,binedges);
    % text2bar(fnames(i).spBothPstSortPstHisto,'',ps(i).lc_bothPstPkUniformity,0.9)

    [fnames(i).spBothPstSortPreFig] = plot_unitWaterfall(alnSpMap(i).pst,binedges,mSort(i).sortPre,1,0);
    plot([mBin mBin],[0 sum(unitgrp)],'r--','LineWidth',2); title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')

    [~,ps(i).pk_distro] = kstest2(nanmean(mousePkPreSum(bvgrp(i).bvInd,:)),nanmean(mousePkPstSum(bvgrp(i).bvInd,:)));
    fnames(i).siBothPkDistroHisto = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.14]); fixRatio(fnames(i).siBothPkDistroHisto);
    [preciup,precidn] = get_CI(mousePkPreSum(bvgrp(i).bvInd,:));
    [pstciup,pstcidn] = get_CI(mousePkPstSum(bvgrp(i).bvInd,:));
    plot_CIs(xcoords,preciup,precidn,fvncols(1,:))
    plot_CIs(xcoords,pstciup,pstcidn,fvncols(2,:))
    plot(xcoords,nanmean(mousePkPreSum(bvgrp(i).bvInd,:)),'Color',fvncols(1,:),'LineWidth',2)
    plot(xcoords,nanmean(mousePkPstSum(bvgrp(i).bvInd,:)),'Color',fvncols(2,:),'LineWidth',2)
    % plot(xcoords,sum(pkMapPre)./sum(pkMapPre,'all'),'Color',vColors2(1,:),'LineWidth',2)
    % plot(xcoords,sum(pkMapPst)./sum(pkMapPst,'all'),'Color',vColors2(2,:),'LineWidth',2)
    plot([0 0],[0 0.3],'k--','LineWidth',2);
    text2bar(fnames(i).siBothPkDistroHisto,'',ps(i).pk_distro);
    set(gca,'Position',[0.11 0.17 0.8 0.80]); xlim([-90 95]); ylim([0 0.3])
    xticks(100*[binedges(1), binedges(round(nBins/2)), binedges(nBins)]);
    set(gca,'FontSize',12,'FontName','Arial')

    if saveFlag
        fsave(fnames(i).spBothPreSortPreFig,[sbase 'sp_' grptag '_pre_sortPre' grp(i).sname])
        fsave(fnames(i).spBothPstSortPstFig,[sbase 'sp_' grptag '_pst_sortPst' grp(i).sname])
        fsave(fnames(i).spBothPstSortPreFig,[sbase 'sp_' grptag '_pst_sortPre' grp(i).sname])
        % fsave(fnames(i).spBothPreSortPreHisto,[sbase 'sp_both_pre_sortPre_hist' grp(i).sname])
        % fsave(fnames(i).spBothPstSortPstHisto,[sbase 'sp_both_pst_sortPst_hist' grp(i).sname])
        fsave(fnames(i).siBothPkDistroHisto,[sbase 'sp_' grptag '_hist' grp(i).sname])
        close all
    end
end

%% Finding Track Relative or Reward Relative cells

for i = 1:2
    spPk1 = binpos(lcDat(:,4));
    spPk2 = binpos(lcDat(:,8));

    dth = 0.3;  % Distance from peak threshold (m)
    drz = r2pos - r1pos;
    cells(i).tr = (spPk2 - spPk1 < dth & spPk2 - spPk1 > -dth)' & siBothID & grp(i).logi; % threshold 30cm
    cells(i).rr = ((spPk2 - spPk1 > drz-dth & spPk2 - spPk1 < drz+dth)' & siBothID  & grp(i).logi) | ((spPk2 - spPk1 < -drz+dth & spPk2 - spPk1 > -drz-dth)' & siBothID  & grp(i).logi);
    cells(i).ir = (siBothID & grp(i).logi) & ~cells(i).tr & ~cells(i).rr;
    xrand = 2*rand(size(lcDat(unitgrp & grp(i).logi,4)))-1;
    yrand = 2*rand(size(lcDat(unitgrp & grp(i).logi,4)))-1;

    pkPosPatchFig = figure; hold on; axis square
    set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.39])
    pkPosPatchFig = fixRatio(pkPosPatchFig);
    patch(100*[0 dth 1.85 1.85 1.85-dth 0 0],100*[0 0 1.85-dth 1.85 1.85 dth 0],'b','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    patch(100*[0 1.25 0.65 0 0],100*[0.6 1.85 1.85 1.2 0.6],'r','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    patch(100*[1.2 1.85 1.85 0.6 1.2],100*[0 0.65 1.25 0 0],'r','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    patch(100*[1.2 1.85 1.85 1.2],100*[0 0 0.65 0],[1 0 1],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    patch(100*[0 0 0.65 0],100*[1.2 1.85 1.85 1.2],[1 0 1],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    patch(100*[dth 0.6 1.85 1.85 dth],100*[0 0 1.25 1.55 0],[1 0 1],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    patch(100*[0 0 1.25 1.55 0],100*[dth 0.6 1.85 1.85 dth],[1 0 1],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(100*spPk1(unitgrp & grp(i).logi)+xrand',100*spPk2(unitgrp & grp(i).logi)+yrand','k.','MarkerSize',10)
    plot([0 100*binpos(end)],[0 100*binpos(end)],'k--')
    % plot([r1pos r1pos]*100,[0 100*binpos(end)],'r--',[0 100*binpos(end)],[r2pos r2pos]*100,'r--')
    xlabel('Absolute Peak Loc. (F)'); xlim([0 100*binpos(end)])
    ylabel('Absolute Peak Loc. (N)'); ylim([0 100*binpos(end)])
    set(gca,'FontSize',16,'FontName','Arial')

    nBoth = sum(siBothID & grp(i).logi);
    nPre  = sum(siFrstID & grp(i).logi);
    nPost = sum(siLastID & grp(i).logi);
    nNot  = sum(not(siFrstID | siLastID) & useCC & grp(i).logi);

    prcts = [sum(cells(i).tr), sum(cells(i).ir), sum(cells(i).rr), nPre - nBoth, nPost - nBoth, nNot] ./ sum(useCC & grp(i).logi);

    cMap = [.75, .75, 1; .75, .5, .75; 1, .75, .75; .35, .35, .35; 1, .25, .25; 1, 1, 1];
    trirrrPie = figure;
    p = piechart(round(prcts*100,2),["TR","IR","RR","F-only","N-only","Neither"]);
    p.LabelStyle = 'namedata';
    colororder(cMap)

    if saveFlag
        fsave(pkPosPatchFig,[sbase 'lc_' grptag '_PeakComp_patch' grp(i).sname])
        fsave(trirrrPie,    [sbase 'lc_PeakComp_pie' grp(i).sname])
    end

end

%% Shuffle novel RZ peak locations and calculate global confidence bands 

for i = 1:2

    activeCells = siBothID & grp(i).logi; %cells(i).rr; % siBothID; rrCells; trCells
    alnPkPst = mod(lcDat(activeCells,8) + rShiftUnit(activeCells), length(binedges)-1)+1;
    alnDistro = histcounts(binedges(alnPkPst),binedges);
    clear deltaField_DistroJit
    for j = 1:250
        rShift = 1 + randi(length(binedges) - 2,sum(activeCells),1);
        jitPks = binedges(mod(alnPkPst + rShift, length(binedges)-1)+1);
        jitPkDistro(:,j) = histcounts(jitPks,binedges);

        % fieldAlignJit = mod(binpos(jitPost)+shiftR2,trackLen)+0.5*dbnsz;
        % fieldDstJit_Align = histcounts(fieldAlignJit,binedges);
        % deltaField_RZJit = fieldAlignJit - fieldAlignR1;
        % circAlignNeg = deltaField_RZJit < -trackLen/2;
        % circAlignPos = deltaField_RZJit > trackLen/2;
        % deltaField_RZJit(circAlignNeg) = -(deltaField_RZJit(circAlignNeg) + trackLen);     % When new field back-shifts
        % deltaField_RZJit(circAlignPos) = -(deltaField_RZJit(circAlignPos) - trackLen);     % When new field forward-shifts
        % deltaField_DistroJit(:,j) = histcounts(deltaField_RZJit,shiftbins);
        % 
        % % uFieldAlignJit = mean([fieldAlignR1; fieldAlignJit],1);
        % % uFieldDistroJit(:,i) = histcounts(uFieldAlignJit - trackLen/2,shiftbins);
        % uFieldDistroJit(:,j) = histcounts(fieldAlignJit-trackLen/2,shiftbins);
    end

    [~,fieldShiftRZJitFig] = get_confband(jitPkDistro',alnDistro,1,binpos*100-92.5,5);
    xlabel('\Delta RZ-aligned Novel - Familiar (cm)');
end
[~,fieldShiftRZJitFig] = get_confband(deltaField_DistroJit',deltaField_Distro,1,shiftbins(1:end-1)*100,dbnsz*100);
xlabel('\Delta RZ-aligned Novel - Familiar (cm)');

% [~,fieldShiftAbsJitFig] = get_confband(uFieldDistroJit',uFieldAlignDistro,1,shiftbins(1:end-1)*100,dbnsz*100);

xcoords = (shiftbins(1:end-1) + 0.5* dbnsz)*100;
d2rzDistroFig = figure; hold on;
set(gcf,'units','normalized','position',[0.3536 0.4231 0.25 0.204])
bar(xcoords,uFieldAlignDistro ./ sum(uFieldAlignDistro));
plot(xcoords,mean(uFieldDistroJit,2) ./ sum(mean(uFieldDistroJit,2)),'k');
ciup = prctile(uFieldDistroJit',95,1) ./ sum(uFieldAlignDistro); % Switch to 97.5 for 2-tail
% cidn = prctile(uFieldDistroJit',2.5,1) ./ sum(uFieldAlignDistro); % For 2-tailed 95% CI
% plot_CIs(xcoords,ciup,cidn,[0 0 0])
plot(xcoords,ciup,'k--')
ylabel('P(RR field peak)'); xlabel('Mean dist. to reward');
legend({'Real','Shuffle Mean','Shuffle 95%'},'location','nw')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(pkPosFig,[sbase 'lc_PeakComp'],'png')
    saveas(pkPosMapFig,[sbase 'lc_PeakMap'],'png')
    saveas(prepstFieldDstFig,[sbase 'lc_Distro'],'png')
    saveas(rzDstFig,[sbase 'lc_Align_Distro'],'png')
    saveas(fieldShiftRZFig,[sbase 'lc_Align_Delta'],'png')
    saveas(fieldShiftRZJitFig,[sbase 'lc_Align_Delta_shuf'],'png')
    fsave(d2rzDistroFig,[sbase 'lc_Align_Dist2RZ_RR'])
end

%% Waterfall by TR or RR

cMap = [0.25 0.15 1; 0.75 0.75 1; 0.25 0.25 0.25];

spTrRrPie = figure;
p = piechart([sum(trCells), sum(rrCells), sum(siBothID)-sum(trCells)-sum(rrCells)],["Track-Rel.","Reward-Rel.","Intermediate"]);
p.LabelStyle = 'namepercent';
colororder(cMap)

[spTRPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(lcMap(trCells,1:length(binpos)),binedges,0,1,0);
plot([r1posInd r1posInd],[0 sum(siBothID)],'k--','LineWidth',2); title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
[spTRPreSortPreHisto,tmpPks] = plot_unitPkHisto(tmpMap,binedges*100,1);
plot([r1pos r1pos]*100,[0 0.11],'k--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.11])
[ps.lc_trPrePkUniformity, stats.lc_trPrePkUniformity] = pkChi2(tmpPks,binedges);
text2bar(spTRPreSortPreHisto,'',ps.lc_trPrePkUniformity,0.9)

[spTRPstSortPreFig,tmpMap] = plot_unitWaterfall(lcMap(trCells,length(binpos)+1:end),binedges,sortPre,1,0);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2); title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')
[spTRPstSortPreHisto,tmpPks] = plot_unitPkHisto(tmpMap,binedges*100,1);
plot([r2pos r2pos]*100,[0 0.11],'r--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.11])
[ps.lc_trPstPkUniformity, stats.lc_trPstPkUniformity] = pkChi2(tmpPks,binedges);
text2bar(spTRPstSortPreHisto,'',ps.lc_trPstPkUniformity,0.9)

[spRRPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(lcMap(rrCells,1:length(binpos)),binedges,0,1,0);
plot([r1posInd r1posInd],[0 sum(siBothID)],'k--','LineWidth',2); title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
[spRRPreSortPreHisto,tmpPks] = plot_unitPkHisto(tmpMap,binedges*100,1);
plot([r1pos r1pos]*100,[0 0.11],'k--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.11])
[ps.lc_rrPrePkUniformity, stats.lc_rrPrePkUniformity] = pkChi2(tmpPks,binedges);
text2bar(spRRPreSortPreHisto,'',ps.lc_rrPrePkUniformity,0.9)

[spRRPstSortPreFig,tmpMap] = plot_unitWaterfall(lcMap(rrCells,length(binpos)+1:end),binedges,sortPre,1,0);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2); title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')
[spRRPstSortPreHisto,tmpPks] = plot_unitPkHisto(tmpMap,binedges*100,1);
plot([r2pos r2pos]*100,[0 0.11],'r--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.11])
[ps.lc_rrPstPkUniformity, stats.lc_rrPstPkUniformity] = pkChi2(tmpPks,binedges);
text2bar(spRRPstSortPreHisto,'',ps.lc_rrPstPkUniformity,0.9)

% [spIRPreSortPreFig,~,sortPre] = plot_unitWaterfall(lcMap(irCells,1:length(binpos)),binedges);
% plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
% 
% spIRPstSortPreFig = plot_unitWaterfall(lcMap(irCells,length(binpos)+1:end),binedges,sortPre);
% plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')

if saveFlag
    fsave(spTrRrPie,[sbase,'sp_trrr_pie'])
    fsave(spTRPreSortPreFig,[sbase 'sp_tr_pre_sortPre'])
    fsave(spTRPstSortPreFig,[sbase 'sp_tr_pst_sortPre'])
    fsave(spRRPreSortPreFig,[sbase 'sp_rr_pre_sortPre'])
    fsave(spRRPstSortPreFig,[sbase 'sp_rr_pst_sortPre'])
    fsave(spTRPreSortPreHisto,[sbase 'sp_tr_pre_sortPre_histo'])
    fsave(spTRPstSortPreHisto,[sbase 'sp_tr_pst_sortPre_histo'])
    fsave(spRRPreSortPreHisto,[sbase 'sp_rr_pre_sortPre_histo'])
    fsave(spRRPstSortPreHisto,[sbase 'sp_rr_pst_sortPre_histo'])
    % saveas(spIRPreSortPreFig,[sbase 'sp_ir_pre_sortPre'],'png')
    % saveas(spIRPstSortPreFig,[sbase 'sp_ir_pst_sortPre'],'png')
end

%% Compare RR and TR parameters pre and post

% FR of TR vs RR cells
[~,ps.fr_PreStnd_trrr,~,stats.fr_PreStnd_trrr] = ttest2(frDat(trCells,1),frDat(rrCells,1));
[~,ps.fr_PreRunn_trrr,~,stats.fr_PreRunn_trrr] = ttest2(frDat(trCells,2),frDat(rrCells,2));
[~,ps.fr_PstStnd_trrr,~,stats.fr_PstStnd_trrr] = ttest2(frDat(trCells,3),frDat(rrCells,3));
[~,ps.fr_PstRunn_trrr,~,stats.fr_PstRunn_trrr] = ttest2(frDat(trCells,4),frDat(rrCells,4));

bardat = [mean(frDat(trCells,1)), mean(frDat(rrCells,1)); mean(frDat(trCells,3)), mean(frDat(rrCells,3));...
    mean(frDat(trCells,2)), mean(frDat(rrCells,2)); mean(frDat(trCells,4)), mean(frDat(rrCells,4))];
sqrr = sqrt(sum(rrCells)); sqtr = sqrt(sum(trCells));
stddat = [std(frDat(trCells,1))./sqtr, std(frDat(rrCells,1))./sqrr; std(frDat(trCells,3))./sqtr, std(frDat(rrCells,3))./sqrr;...
    std(frDat(trCells,2))./sqtr, std(frDat(rrCells,2))./sqrr; std(frDat(trCells,4))./sqtr, std(frDat(rrCells,4))./sqrr];

frSTRRRFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
b1 = bar(bardat,'FaceColor','flat');
b1(1).CData = [0 0 1];
b1(2).CData = [1 0 0];
errorbar([0.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15],reshape(bardat',8,1),reshape(stddat',8,1),'k.')
xlim([0.5 4.5]); xticks(1:4); xticklabels({'Stand Fam.', 'Stand Nov.', 'Run Fam.', 'Run Nov.'})
ylabel('Firing Rate (Hz)'); legend({'TR','RR'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

% Spatial info
[~,ps.si_pre_trrr,~,stats.si_pre_trrr] = ttest2(lcDat(trCells,2),lcDat(rrCells,2));
[~,ps.si_pst_trrr,~,stats.si_pst_trrr] = ttest2(lcDat(trCells,6),lcDat(rrCells,6));

siTRRRPreFig = plotMiniBar(lcDat(trCells,2),lcDat(rrCells,2),[0 0 1; 1 0 0]);
xticklabels({'TR','RR'}); ylim([0 5]); text2bar(siTRRRPreFig,'Spatial Info.',ps.si_pre_trrr);
siTRRRPstFig = plotMiniBar(lcDat(trCells,6),lcDat(rrCells,6),[0 0 1; 1 0 0]);
xticklabels({'TR','RR'}); ylim([0 5]); text2bar(siTRRRPstFig,'Spatial Info.',ps.si_pst_trrr);

% SI Delta histogram
[~,ps.si_dlt_tr,~,stats.si_dlt_tr] = ttest(lcDat(trCells & siBothID,2) - lcDat(trCells & siBothID,6));
[~,ps.si_dlt_rr,~,stats.si_dlt_rr] = ttest(lcDat(rrCells & siBothID,2) - lcDat(rrCells & siBothID,6));
binedges = -2:0.2:2;
siTRRRDltFig = plotDeltaHisto2(lcDat(trCells,2), lcDat(trCells,6),...
    lcDat(rrCells,2), lcDat(rrCells,6), [ps.si_dlt_tr, ps.si_dlt_rr], binedges, [0 0 1; 1 0 0]);
xlim(binedges([1,end])); xlabel('\Delta Spatial Info (N - F)')

% Theta MRL and angle
[~,ps.thA_pre_trrr,~,stats.thA_pre_trrr] = ttest2(thDat(trCells,3),thDat(rrCells,3));
[~,ps.thA_pst_trrr,~,stats.thA_pst_trrr] = ttest2(thDat(trCells,6),thDat(rrCells,6));
[~,ps.thM_pre_trrr,~,stats.thM_pre_trrr] = ttest2(thDat(trCells,2),thDat(rrCells,2));
[~,ps.thM_pst_trrr,~,stats.thM_pst_trrr] = ttest2(thDat(trCells,5),thDat(rrCells,5));

thATRRRPreFig = plotMiniBar(thDat(trCells,3),thDat(rrCells,3),[0 0 1; 1 0 0]);
xticklabels({'TR','RR'}); ylim([-180 180]); text2bar(thATRRRPreFig,'Theta angle',ps.thA_pre_trrr);
thATRRRPstFig = plotMiniBar(thDat(trCells,6),thDat(rrCells,6),[0 0 1; 1 0 0]);
xticklabels({'TR','RR'}); ylim([-180 180]); text2bar(thATRRRPstFig,'Theta angle',ps.thA_pst_trrr);
% thMTRRRPreFig = plotMiniBar(thDat(trCells,2),thDat(rrCells,2),[0 0 1; 1 0 0]);
% xticklabels({'TR','RR'}); text2bar(thMTRRRPreFig,'Theta MRL',ps.thM_pre_trrr);
% thMTRRRPstFig = plotMiniBar(thDat(trCells,4),thDat(rrCells,4),[0 0 1; 1 0 0]);
% xticklabels({'TR','RR'}); text2bar(thMTRRRPstFig,'Theta MRL',ps.thM_pst_trrr);

% Theta Delta histograms
[~,ps.thA_dlt_tr,~,stats.thA_dlt_tr] = ttest(thDat(trCells & thBothID,3) - thDat(trCells & thBothID,6));
[~,ps.thA_dlt_rr,~,stats.thA_dlt_rr] = ttest(thDat(rrCells & thBothID,3) - thDat(rrCells & thBothID,6));
binedges = rad2deg(-pi/2:pi/36:pi/2);
thATRRRDltFig = plotDeltaHisto2(thDat(trCells & thBothID,3), thDat(trCells & thBothID,6),...
    thDat(rrCells & thBothID,3), thDat(rrCells & thBothID,6), [ps.thA_dlt_tr, ps.thA_dlt_rr], binedges, [0 0 1; 1 0 0]);
xlim([-50 50]); xlabel('\Delta Theta Angle (N - F)')

% Velocity coding
[~,ps.vlB_pre_trrr,~,stats.vlB_pre_trrr] = ttest2(vlDat(trCells,2),vlDat(rrCells,2));
[~,ps.vlB_pst_trrr,~,stats.vlB_pst_trrr] = ttest2(vlDat(trCells,5),vlDat(rrCells,5));

vlTRRRPreFig = plotMiniBar(vlDat(trCells,2),vlDat(rrCells,2),[0 0 1; 1 0 0]);
xticklabels({'TR','RR'}); ylim([-0.5 1]); text2bar(vlTRRRPreFig,'Velocity Slope',ps.vlB_pre_trrr);
vlTRRRPstFig = plotMiniBar(vlDat(trCells,5),vlDat(rrCells,5),[0 0 1; 1 0 0]);
xticklabels({'TR','RR'}); ylim([-0.5 1]); text2bar(vlTRRRPstFig,'Velocity Slope',ps.vlB_pst_trrr);

% Vl Delta histogram
[~,ps.vlB_dlt_tr,~,stats.vlB_dlt_tr] = ttest(vlDat(trCells & vlBothID,2) - vlDat(trCells & vlBothID,5));
[~,ps.vlB_dlt_rr,~,stats.vlB_dlt_rr] = ttest(vlDat(rrCells & vlBothID,2) - vlDat(rrCells & vlBothID,5));
binedges = -.5:0.05:.5;
vlBTRRRDltFig = plotDeltaHisto2(vlDat(trCells & vlBothID,2), vlDat(trCells & vlBothID,5),...
    vlDat(rrCells & vlBothID,2), vlDat(rrCells & vlBothID,5), [ps.vlB_dlt_tr, ps.vlB_dlt_rr], binedges, [0 0 1; 1 0 0]);
xlim(binedges([1,end])); xlabel('\Delta Velocity Slope (N - F)')

if saveFlag
    fsave(siTRRRPreFig,[sbase 'si_prebar_trrr'])
    fsave(siTRRRPstFig,[sbase 'si_pstbar_trrr'])
    fsave(siTRRRDltFig,[sbase 'si_dlta_trrr'])
    fsave(thATRRRPreFig,[sbase 'thA_prebar_trrr'])
    fsave(thATRRRPstFig,[sbase 'thA_pstbar_trrr'])
    fsave(thATRRRDltFig,[sbase 'thA_dlta_trrr'])
    % fsave(thMTRRRPreFig,[sbase 'thM_prebar_trrr'])
    % fsave(thMTRRRPstFig,[sbase 'thM_pstbar_trrr'])
    fsave(vlTRRRPreFig,[sbase 'vl_prebar_trrr'])
    fsave(vlTRRRPstFig,[sbase 'vl_pstbar_trrr'])
    fsave(vlBTRRRDltFig,[sbase 'vlB_dlta_trrr'])
    fsave(frSTRRRFig,[sbase 'fr_bar_trrr'])
end

%% Population Vector analysis

for i = 1:2
    activeUnits = siBothID & grp(i).logi;     % siBothID (allSICells); trCells; rrCells;
    nBins = length(binpos);

    % posNormPre = lcMap(activeUnits,1:nBins) ./ max(lcMap(activeUnits,1:nBins),[],2);  % Using TR rate map
    % posNormPst = lcMap(activeUnits,nBins+1:end) ./ max(lcMap(activeUnits,nBins+1:end),[],2);  % Using TR rate map
    posNormPre = alnSpMap(i).pre ./ max(alnSpMap(i).pre,[],2);  % RR rate maps, make sure to build them with siBothID in the block above
    posNormPst = alnSpMap(i).pst ./ max(alnSpMap(i).pst,[],2);  % RR rate maps

    pvPrePst = corr(posNormPst,posNormPre);
    fnames(i).pvPrePstF = plot_pvcorr(pvPrePst,[-90 0 90]); ylabel(''); xlabel('');
    if saveFlag
        fsave(fnames(i).pvPrePstF,[sbase 'lc_pv_corr_allSICells_prepost_align' grp(i).sname])
    end
end

%% LME version of diagonal / off diagonal comparison
idMat = logical(eye(size(pvPrePst)));
dgMat = logical(spdiags([1 1],[-round(nBins/2) round(nBins/2)],nBins,nBins));

for i = 1:nMice
    posNormPre = alnMouseMap(i).pre ./ max(alnMouseMap(i).pre,[],2);
    posNormPst = alnMouseMap(i).pst ./ max(alnMouseMap(i).pst,[],2);
    try
        pvPrePst = corr(posNormPst,posNormPre);
    catch
        pvPrePst = NaN(nBins,nBins);
    end
    uPVCorrID(i) = mean(pvPrePst(idMat),'all');
    uPVCorrDG(i) = mean(pvPrePst(dgMat),'all');

    try
        pvOddEvn(i).pre = corr(pvStr(i).preOddsub, pvStr(i).preEvnsub);
        pvOddEvn(i).pst = corr(pvStr(i).pstOddsub, pvStr(i).pstEvnsub);
    catch
        pvOddEvn(i).pre = NaN(nBins,nBins);
        pvOddEvn(i).pst = NaN(nBins,nBins);
    end

    uPVCPreOddEvnID(i) = mean(pvOddEvn(i).pre(idMat),'all','omitnan');
    uPVCPreOddEvnDG(i) = mean(pvOddEvn(i).pre(dgMat),'all','omitnan');
    uPVCPstOddEvnID(i) = mean(pvOddEvn(i).pst(idMat),'all','omitnan');
    uPVCPstOddEvnDG(i) = mean(pvOddEvn(i).pst(dgMat),'all','omitnan');
end

% for j = 1:2
%     [~,ps(j).lc_pv_PrePst_idVdg,~,stats(j).lc_pv_PrePst_idVdg] = ttest(uPVCorrID(bvgrp(j).bvInd),uPVCorrDG(bvgrp(j).bvInd));
%     [~,ps(j).lc_pv_PreOddEvn_idVdg,~,stats(j).lc_pv_PreOddEvn_idVdg] = ttest(uPVCPreOddEvnID(bvgrp(j).bvInd),uPVCPreOddEvnDG(bvgrp(j).bvInd));
%     [~,ps(j).lc_pv_PstOddEvn_idVdg,~,stats(j).lc_pv_PstOddEvn_idVdg] = ttest(uPVCPstOddEvnID(bvgrp(j).bvInd),uPVCPstOddEvnDG(bvgrp(j).bvInd));
% 
%     fnames(j).pvPrePstCompF    = plot_PVCorrComp(uPVCorrID(bvgrp(j).bvInd),uPVCorrDG(bvgrp(j).bvInd), ps(j).lc_pv_PrePst_idVdg,lnlcols(j,:)); xticklabels({'RR','TR'})
%     fnames(j).pvPreOddEvnCompF = plot_PVCorrComp(uPVCPreOddEvnID(bvgrp(j).bvInd),uPVCPreOddEvnDG(bvgrp(j).bvInd),ps(j).lc_pv_PreOddEvn_idVdg,lnlcols(j,:)./2);
%     fnames(j).pvPstOddEvnCompF = plot_PVCorrComp(uPVCPstOddEvnID(bvgrp(j).bvInd),uPVCPstOddEvnDG(bvgrp(j).bvInd),ps(j).lc_pv_PstOddEvn_idVdg,lnlcols(j,:));
% 
%     if saveFlag
%         fsave(fnames(j).pvPrePstCompF,[sbase 'lc_pv_comp_allSICells_prepost_align' grp(j).sname])
%         fsave(fnames(j).pvPreOddEvnCompF,[sbase 'lc_pv_comp_allSICells_preoddeven' grp(j).sname])
%         fsave(fnames(j).pvPstOddEvnCompF,[sbase 'lc_pv_comp_allSICells_pstoddeven' grp(j).sname])
%         close all
%     end
% end

pvcVarNames = {'pvc','onIDLine','grp','mouse'};
[uPV_FN_lme_table] = get_lmetable([uPVCorrDG'; uPVCorrID'],bvgrp,mID,pvcVarNames);
uPV_FN_lme = fitlme(uPV_FN_lme_table,'pvc ~ onIDLine * grp + (1|mouse)');
pv_FN_lmeF = plot_2wayLME(uPVCorrID,uPVCorrDG,bvgrp,lnlcols);
xticklabels({'RR','TR','RR','TR'})

[uPV_FOE_lme_table] = get_lmetable([uPVCPreOddEvnDG'; uPVCPreOddEvnID'],bvgrp,mID,pvcVarNames);
uPV_FOE_lme = fitlme(uPV_FOE_lme_table,'pvc ~ onIDLine * grp + (1|mouse)');
pv_FOE_lmeF = plot_2wayLME(uPVCPreOddEvnID,uPVCPreOddEvnDG,bvgrp,lnlcols);
xticklabels({'ID','Off-ID','ID','Off-ID'})

[uPV_NOE_lme_table] = get_lmetable([uPVCPstOddEvnDG'; uPVCPstOddEvnID'],bvgrp,mID,pvcVarNames);
uPV_NOE_lme = fitlme(uPV_NOE_lme_table,'pvc ~ onIDLine * grp + (1|mouse)');
pv_NOE_lmeF = plot_2wayLME(uPVCPstOddEvnID,uPVCPstOddEvnDG,bvgrp,lnlcols);
xticklabels({'ID','Off-ID','ID','Off-ID'})

if saveFlag
    fsave(pv_FN_lmeF,[sbase 'lc_pv_comp_allSICells_prepost_align_lme'])
    fsave(pv_FOE_lmeF,[sbase 'lc_pv_comp_allSICells_preoddeven_lme'])
    fsave(pv_NOE_lmeF,[sbase 'lc_pv_comp_allSICells_pstoddeven_lme'])
end

%% PV across laps
PVxLap_FvF = [];
PVxLap_NvN = [];

for i = 1:nMice
    try
        PVxLap_FvF(:,i) = pvStr(i).preBlockPVsub(1:50,1);
    catch
        PVxLap_FvF(:,i) = NaN;
    end
    try
        PVxLap_NvN(:,i) = pvStr(i).pstBlockPVsub(1:50,1);
    catch
        PVxLap_NvN(:,i) = NaN;
    end
end

mdl_PVxLap_FvF_ln = get_linfit(xsF,mean(PVxLap_FvF(:,bvgrp(1).bvInd),2,'omitmissing'));
mdl_PVxLap_FvF_nl = get_linfit(xsF,mean(PVxLap_FvF(:,bvgrp(2).bvInd),2,'omitmissing'));
mdl_PVxLap_NvN_ln = get_linfit(xsF,mean(PVxLap_NvN(:,bvgrp(1).bvInd),2,'omitmissing'));
mdl_PVxLap_NvN_nl = get_linfit(xsF,mean(PVxLap_NvN(:,bvgrp(2).bvInd),2,'omitmissing'));

[ciup_FvF_ln, cidn_FvF_ln] = get_CI(PVxLap_FvF(:,bvgrp(1).bvInd)');
[ciup_FvF_nl, cidn_FvF_nl] = get_CI(PVxLap_FvF(:,bvgrp(2).bvInd)');
[ciup_NvN_ln, cidn_NvN_ln] = get_CI(PVxLap_NvN(:,bvgrp(1).bvInd)');
[ciup_NvN_nl, cidn_NvN_nl] = get_CI(PVxLap_NvN(:,bvgrp(2).bvInd)');

pvXt_lnF = figure; hold on;
set(gcf,'Units','normalized','Position',[1.2 0.4 0.6315 0.1703])
plot_CIs(xsF, ciup_FvF_ln, cidn_FvF_ln, lnlcols(1,:)/2);
plot_CIs(xsN, ciup_NvN_ln, cidn_NvN_ln, lnlcols(1,:));
plot(xsF, mean(PVxLap_FvF(:,bvgrp(1).bvInd),2,'omitmissing'),'Color',lnlcols(1,:)/2,'LineWidth',2);
plot(xsN, mean(PVxLap_NvN(:,bvgrp(1).bvInd),2,'omitmissing'),'Color',lnlcols(1,:),'LineWidth',2);
plot(xsF,mdl_PVxLap_FvF_ln.ypred, 'color', lnlcols(1,:)/2, 'LineWidth',2)
plot(xsN,mdl_PVxLap_NvN_ln.ypred, 'color', lnlcols(1,:),   'LineWidth',2)
ylim([0 1]); xlabel('Laps')
text2bar(pvXt_lnF,"",    mdl_PVxLap_FvF_ln.p, 0.8, 0.1, lnlcols(1,:)/2);
text2bar(pvXt_lnF,"PVC", mdl_PVxLap_NvN_ln.p, 0.3, 0.4, lnlcols(1,:));
set(gca,'FontSize',16,'FontName','Arial')

pvXt_nlF = figure; hold on;
set(gcf,'Units','normalized','Position',[1.2 0.4 0.6315 0.1703])
plot_CIs(xsF, ciup_FvF_nl, cidn_FvF_nl, lnlcols(2,:)/2);
plot_CIs(xsN, ciup_NvN_nl, cidn_NvN_nl, lnlcols(2,:));
plot(xsF, mean(PVxLap_FvF(:,bvgrp(2).bvInd),2,'omitmissing'),'Color',lnlcols(2,:)/2,'LineWidth',2);
plot(xsN, mean(PVxLap_NvN(:,bvgrp(2).bvInd),2,'omitmissing'),'Color',lnlcols(2,:),'LineWidth',2);
plot(xsF,mdl_PVxLap_FvF_nl.ypred, 'color', lnlcols(2,:)/2, 'LineWidth',2)
plot(xsN,mdl_PVxLap_NvN_nl.ypred, 'color', lnlcols(2,:),   'LineWidth',2)
ylim([0 1]); xlabel('Laps')
text2bar(pvXt_nlF,"",    mdl_PVxLap_FvF_nl.p, 0.8, 0.1, lnlcols(2,:)/2);
text2bar(pvXt_nlF,"PVC", mdl_PVxLap_NvN_nl.p, 0.3, 0.4, lnlcols(2,:));
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    fsave(pvXt_lnF, [sbase 'pv_corrXlap_ln'])
    fsave(pvXt_nlF, [sbase 'pv_corrXlap_nl'])
end

%% Correlate PV to LDI by lap
lapcutoff = 50;
for i = 1:nMice
    nTrialsPre(i) = size(pvStr(i).preBlockPVsub,1);
    nTrialsPst(i) = size(pvStr(i).pstBlockPVsub,1);
end
maxTr = [max(nTrialsPre) max(nTrialsPst)];
pvXlapPrePre = NaN(maxTr(1),nMice);
pvXlapPrePst = NaN(maxTr(1),nMice);
pvXlapPstPre = NaN(maxTr(2),nMice);
pvXlapPstPst = NaN(maxTr(2),nMice);

for i = 1:nMice
    pvXlapPrePre(1:nTrialsPre(i),i) = pvStr(i).preBlockPVsub(:,1);
    pvXlapPrePst(1:nTrialsPre(i),i) = pvStr(i).preBlockPVsub(:,2);
    pvXlapPstPre(1:nTrialsPst(i),i) = pvStr(i).pstBlockPVsub(:,1);
    pvXlapPstPst(1:nTrialsPst(i),i) = pvStr(i).pstBlockPVsub(:,2);
end

% cmapcool = cool(nMice);
% pvXlapFig = figure; hold on;
% set(gcf,'units','normalized','position',[0.4 0.35 0.35 0.35])
% for i = 1:nMice
%     plot(1:maxTr(1),pvXlapPrePre(:,i),'Color',cmapcool(i,:))
% end
% plot(1:maxTr(1),nanmean(pvXlapPrePre'),'k','LineWidth',2)
% for i = 1:nMice
%     plot(101:100+maxTr(2),pvXlapPstPre(:,i),'Color',cmapcool(i,:))
% end
% ylim([-0.25 1]); ylabel('PV Corr. on-diagonal'); xlabel('Lap #')
% plot(101:100+maxTr(2),nanmean(pvXlapPstPre'),'k','LineWidth',2)
% xticks(00:20:200)
% xticklabels([00:20:80,0:20:100])
% legend('Mice')
% set(gca,'FontSize',16,'FontName','Arial')

pvXlapPrePre = pvXlapPrePre(1:lapcutoff,:);
pvXlapPstPst = pvXlapPstPst(1:lapcutoff,:);

pvXlapFig = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.30])
[preCIup,preCIdn] = get_CI(pvXlapPrePre'); 
[pstCIup,pstCIdn] = get_CI(pvXlapPstPst'); 
plot_CIs(1:lapcutoff,preCIup,preCIdn,fvncols(1,:)); % 1:maxTr(1)
plot_CIs(1:lapcutoff,pstCIup,pstCIdn,fvncols(2,:)); % 1:maxTr(2)
plot(1:lapcutoff,nanmean(pvXlapPrePre'),'Color',fvncols(1,:),'LineWidth',2)
plot(1:lapcutoff,nanmean(pvXlapPstPst'),'Color',fvncols(2,:),'LineWidth',2)
ylim([0 1]); ylabel('PV Corr.'); xlabel('Lap #'); xlim([0 lapcutoff]);
legend({'F-->F','N-->N'},'location','nw')
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    fsave(pvXlapFig,[sbase 'lc_pv_corr_laps'])
end

%% Velocity Corr across time

lapcutoff = 65;
for i = 1:nMice
    nTrialsPre(i) = size(bvDat(i).vCorPre,1);
    nTrialsPst(i) = size(bvDat(i).vCorPst,1);
end
maxTr = [max(nTrialsPre) max(nTrialsPst)];
vcXlapPre = NaN(maxTr(1),nMice);
vcXlapPst = NaN(maxTr(1),nMice);

for i = 1:nMice
    vcXlapPre(1:nTrialsPre(i),i) = bvDat(i).vCorPre(:,1);
    vcXlapPst(1:nTrialsPst(i),i) = bvDat(i).vCorPst(:,1);
end

vcXlapPre = vcXlapPre(1:lapcutoff,:);
vcXlapPst = vcXlapPst(1:lapcutoff,:);

vcXlapFig = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.30])
[preCIup,preCIdn] = get_CI(vcXlapPre'); 
[pstCIup,pstCIdn] = get_CI(vcXlapPst'); 
plot_CIs(1:lapcutoff,preCIup,preCIdn,fvncols(1,:))
plot_CIs(1:lapcutoff,pstCIup,pstCIdn,fvncols(2,:))
plot(1:lapcutoff,nanmean(vcXlapPre'),'Color',fvncols(1,:),'LineWidth',2)
plot(1:lapcutoff,nanmean(vcXlapPst'),'Color',fvncols(2,:),'LineWidth',2)
ylim([-0.05 1]); ylabel('Velocity Corr.'); xlabel('Lap #'); xlim([0 lapcutoff]);
legend({'F-->F','N-->N'},'location','nw')
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    fsave(vcXlapFig,[sbase 'bhv_v_corr_laps'])
end

%% Lick DI across time

for i = 1:nMice
    nTrialsPre(i) = size(bvDat(i).preLckDI,1);
    nTrialsPst(i) = size(bvDat(i).pstLckDI,1);
end
maxTr = [max(nTrialsPre) max(nTrialsPst)];
lDIXlapPre = NaN(maxTr(1),nMice);
lDIXlapPst = NaN(maxTr(1),nMice);

for i = 1:nMice
    lDIXlapPre(1:nTrialsPre(i),i) = bvDat(i).preLckDI(:,1);
    lDIXlapPst(1:nTrialsPst(i),i) = bvDat(i).pstLckDI(:,1);
end

lDIXlapPre = lDIXlapPre(1:lapcutoff,:);
lDIXlapPst = lDIXlapPst(1:lapcutoff,:);

lDIXlapFig = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.30])
[preCIup,preCIdn] = get_CI(lDIXlapPre'); 
[pstCIup,pstCIdn] = get_CI(lDIXlapPst'); 
plot_CIs(1:lapcutoff,preCIup,preCIdn,fvncols(1,:))
plot_CIs(1:lapcutoff,pstCIup,pstCIdn,fvncols(2,:))
plot(1:lapcutoff,nanmean(lDIXlapPre'),'Color',fvncols(1,:),'LineWidth',2)
plot(1:lapcutoff,nanmean(lDIXlapPst'),'Color',fvncols(2,:),'LineWidth',2)
ylim([-0.05 1]); ylabel('Lick DI'); xlabel('Lap #'); xlim([0 lapcutoff]);
legend({'Familiar','Novel'},'location','nw')
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    fsave(lDIXlapFig,[sbase 'bhv_lickDI_laps'])
end

%% Velocity vs PV Corr analysis

vcXpvcPrePre = [];
vcXpvcPrePst = [];
vcXpvcPstPre = [];
vcXpvcPstPst = [];

for i = 1:nMice
    % vcXpvcPrePre = [vcXpvcPrePre; vcXlapPre(:,i), pvXlapPrePre(2:end,i)];   % [Vcorr, PVC] concatenated for all mice
    % vcXpvcPrePst = [vcXpvcPrePst; vcXlapPre(:,i), pvXlapPrePst(2:end,i)];
    vcXpvcPrePre = [vcXpvcPrePre; vcXlapPre(:,i), pvXlapPrePre(:,i)];   % [Vcorr, PVC] concatenated for all mice
    % vcXpvcPrePst = [vcXpvcPrePst; vcXlapPre(:,i), pvXlapPrePst(:,i)];
    % vcXpvcPstPre = [vcXpvcPstPre; vcXlapPst(:,i), pvXlapPstPre(:,i)];
    vcXpvcPstPst = [vcXpvcPstPst; vcXlapPst(:,i), pvXlapPstPst(:,i)];
end
preNans = isnan(vcXpvcPrePre); preNans = sum(preNans,2) > 0;
pstNans = isnan(vcXpvcPstPst); pstNans = sum(pstNans,2) > 0;
vcXpvcPrePre(preNans,:) = [];
vcXpvcPstPst(pstNans,:) = [];

mdlVPreXPVCprepre = get_linfit(vcXpvcPrePre(:,1),vcXpvcPrePre(:,2));
mdlVPreXPVCpstpst = get_linfit(vcXpvcPstPst(:,1),vcXpvcPstPst(:,2));

vcXpvcF = figure; hold on
plot(vcXpvcPrePre(:,1),vcXpvcPrePre(:,2),'o','Color',fvncols(1,:))
plot(vcXpvcPstPst(:,1),vcXpvcPstPst(:,2),'o','Color',fvncols(2,:))
plot(vcXpvcPrePre(:,1),mdlVPreXPVCprepre.ypred,'Color',fvncols(1,:),'LineWidth',2)
plot(vcXpvcPstPst(:,1),mdlVPreXPVCpstpst.ypred,'Color',fvncols(2,:),'LineWidth',2)
xlim([-1 1]); xlabel('Velocity Correlation');
ylim([0 1]); ylabel('PV Correlation');
legend({'F-F Vel -> F-F PVC','N-N Vel -> N-N PVC'},'location','ne')
set(gca,'FontSize',16,'FontName','Arial')
xlims = xlim;
ylims = ylim;
text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims),  ['R = ' num2str(mdlVPreXPVCprepre.r, 3)], 'Color', fvncols(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlVPreXPVCprepre.p, 3)], 'Color', fvncols(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.25*diff(ylims), ['R = ' num2str(mdlVPreXPVCpstpst.r, 3)], 'Color', fvncols(2,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.3*diff(ylims),  ['p = ' num2str(mdlVPreXPVCpstpst.p, 3)], 'Color', fvncols(2,:), 'FontSize', 12)

if saveFlag
    fsave(vcXpvcF,[sbase '_vCorrXPVC'])
end

%% SPWR Data
% rpDat: 2&4 = sig.; 1&3 = participation; 
% rpRat: 1&4 = rate by mouse; 2&5 = mean duration; 3&6 = %age ripples > 100ms
excludemice = [32 79];

for i = 1:2
    % Quantify proportions of spatially modulated units
    [fnames(i).swrPie, fnames(i).swrPrcts] = prepostPie(swrFrstID(grp(i).inds),swrLastID(grp(i).inds),useCC(grp(i).inds));
    frstGrpAdd = [mID(bvgrp(i).bvInd); recID(swrFrstID & grp(i).logi,1)]; % Need to add all mice then subtract 1 from groupCounts otherwise it is blind to missing mice
    lastGrpAdd = [mID(bvgrp(i).bvInd); recID(swrLastID & grp(i).logi,1)];
    swrBothCounts = [groupcounts(frstGrpAdd)-1, groupcounts(lastGrpAdd)-1];
    swrBothRatio(i).ratio = swrBothCounts ./ repmat(groupcounts([mID(bvgrp(i).bvInd); recID(useCC & grp(i).logi,1)])-1,[1,2]);
    swrBothRatio(i).ratio(isnan(swrBothRatio(i).ratio(:,1)),:) = [];
    [~,ps(i).swr_ModCt_both,~,stats(i).swr_ModCt_both] = ttest(swrBothCounts(:,1),swrBothCounts(:,2));
    fnames(i).swrModCtFig = plot_barXmouse(swrBothRatio(i).ratio);
    text2bar(fnames(i).swrModCtFig,'P(Sig. SWR Mod)',ps(i).swr_ModCt_both);

    clear useMice
    for j = 1:length(bvgrp(i).mID)
        if sum(excludemice == bvgrp(i).mID(j)) > 0 % Remove mice with only CA1 data
            useMice(j) = 0;
        else
            useMice(j) = 1;
        end
    end
    useMice = logical(useMice);
    tmpRpDat = rpRat(bvgrp(i).bvInd,:);

    [~,ps(i).swr_SWRRate,~,stats(i).swr_SWRRate]   = ttest(tmpRpDat(useMice,1),tmpRpDat(useMice,4));
    [~,ps(i).swr_SWRDur,~,stats(i).swr_SWRDur]     = ttest(tmpRpDat(useMice,2),tmpRpDat(useMice,5));
    [~,ps(i).swr_SWRLng,~,stats(i).swr_SWRLng]     = ttest(tmpRpDat(useMice,3),tmpRpDat(useMice,6));
    [~,ps(i).swr_PPP_both,~,stats(i).swr_PPP_both] = ttest(rpDat(swrBothID & grp(i).logi,1),rpDat(swrBothID & grp(i).logi,3));
    [~,ps(i).swr_PPP_eith,~,stats(i).swr_PPP_eith] = ttest2(rpDat(swrFrstID & grp(i).logi,1),rpDat(swrLastID & grp(i).logi,3));

    % Test ripple parameters - rate, duration, and %age >100ms
    fnames(i).rpRatFig = plot_barXmouse(tmpRpDat(useMice,[1 4]));
    text2bar(fnames(i).rpRatFig,'SPW-R Rate (Hz)',ps(i).swr_SWRRate);
    fnames(i).rpDurFig = plot_barXmouse(tmpRpDat(useMice,[2 5]));
    text2bar(fnames(i).rpDurFig,'SPW-R Duration (ms)',ps(i).swr_SWRDur);
    fnames(i).rpLngFig = plot_barXmouse(tmpRpDat(useMice,[3 6]));
    text2bar(fnames(i).rpLngFig,'P(SPW-R > 100ms)',ps(i).swr_SWRLng);

    % Units modulated in both phases
    fnames(i).swPartcpBothFig = plotBar2(rpDat(swrBothID,1),rpDat(swrBothID,3));
    text2bar(fnames(i).swPartcpBothFig,'P(SPW-R Participation)',ps(i).swr_PPP_both);

    % For units modulated in either task phase
    fnames(i).swPartcpEithFig = plotBar2(rpDat(swrFrstID & ~swrLastID,1),rpDat(swrLastID & ~swrFrstID,3));
    text2bar(fnames(i).swPartcpEithFig,'P(SPW-R Participation)',ps(i).swr_PPP_eith);

    if saveFlag
        fsave(fnames(i).swrPie,[sbase 'swr_Mod_pie'])
        fsave(fnames(i).rpRatFig,[sbase 'swr_Rate_bar'])
        fsave(fnames(i).rpDurFig,[sbase 'swr_duration_bar'])
        fsave(fnames(i).rpLngFig,[sbase 'swr_longrip_bar'])
        fsave(fnames(i).swPartcpBothFig,[sbase 'swr_Partcp_bar'])
        fsave(fnames(i).swrModCtFig,[sbase 'swr_ModCt_bar'])
        fsave(fnames(i).swPartcpEithFig,[sbase 'swr_Partcp_eith'])
        close all
    end
end

% [~,ps(i).swr_PPPre_trrr,~,stats(i).swr_PPPre_trrr] = ttest2(rpDat(trCells,1),rpDat(rrCells,1));
% [~,ps(i).swr_PPPst_trrr,~,stats(i).swr_PPPst_trrr] = ttest2(rpDat(trCells,3),rpDat(rrCells,3));

%% SPWR lme analysis

cleanRpRat = rpRat;
cleanMID   = mID;
cleanbhvID = bhvID;
for i = 1:length(excludemice)
    rmInds(i) = find(mID == excludemice(i));
end
cleanRpRat(rmInds,:) = []; % Remove those rows
cleanMID(rmInds,:)   = [];
cleanbhvID(rmInds,:) = [];
mLern = [101 99 97 77 73 35 6];
mNonl = [100 91 87 82 80 29 20 12];

lnInd = [];
nlInd = [];

for i = 1:size(cleanbhvID,1)
    if ~isempty(find(mLern == cleanbhvID(i,1), 1))
        lnInd = [lnInd; find(cleanbhvID(:,1) == cleanbhvID(i,1))];
    else
        nlInd = [nlInd; find(cleanbhvID(:,1) == cleanbhvID(i,1))];
    end
end

cleanbvgrp(1).bvInd = lnInd;
cleanbvgrp(2).bvInd = nlInd;
cleanbvgrp(1).grpname = 'learn';
cleanbvgrp(2).grpname = 'nolearn';
cleanbvgrp(1).mID = cleanbhvID(cleanbvgrp(1).bvInd);
cleanbvgrp(2).mID = cleanbhvID(cleanbvgrp(2).bvInd);
cleanbvgrp(1).n = length(cleanbvgrp(1).bvInd);
cleanbvgrp(2).n = length(cleanbvgrp(2).bvInd);

swrVarNames = {'ripRate','isNovel','grp','mouse'};
[swrRat_FN_lme_table] = get_lmetable([cleanRpRat(:,1); cleanRpRat(:,4)],cleanbvgrp,cleanMID,swrVarNames);
swrRat_FN_lme = fitlme(swrRat_FN_lme_table,'ripRate ~ isNovel * grp + (1|mouse)');
swrRat_FN_lmeF = plot_2wayLME(cleanRpRat(:,1)',cleanRpRat(:,4)',cleanbvgrp,lnlcols);
xticklabels({'F','N','F','N'}); ylim([0 2]); ylabel('SPW-R Rate (Hz)')

swrVarNames = {'ripDur','isNovel','grp','mouse'};
[swrDur_FN_lme_table] = get_lmetable([cleanRpRat(:,2); cleanRpRat(:,5)],cleanbvgrp,cleanMID,swrVarNames);
swrDur_FN_lme = fitlme(swrDur_FN_lme_table,'ripDur ~ isNovel * grp + (1|mouse)');
swrDur_FN_lmeF = plot_2wayLME(cleanRpRat(:,2)',cleanRpRat(:,5)',cleanbvgrp,lnlcols);
xticklabels({'F','N','F','N'}); ylim([0 130]); ylabel('SPW-R Dur. (ms)')

swrVarNames = {'ripLng','isNovel','grp','mouse'};
[swrLng_FN_lme_table] = get_lmetable([cleanRpRat(:,3); cleanRpRat(:,6)],cleanbvgrp,cleanMID,swrVarNames);
swrLng_FN_lme = fitlme(swrLng_FN_lme_table,'ripLng ~ isNovel * grp + (1|mouse)');
swrLng_FN_lmeF = plot_2wayLME(cleanRpRat(:,3)',cleanRpRat(:,6)',cleanbvgrp,lnlcols);
xticklabels({'F','N','F','N'}); ylim([0 1]); ylabel('P(SPW-R > 100ms)')

ripPtp = zeros(length(cleanMID),2);
ripPtp(cleanbvgrp(1).bvInd,:) = swrBothRatio(1).ratio;
ripPtp(cleanbvgrp(2).bvInd,:) = swrBothRatio(2).ratio;
swrVarNames = {'ripPtp','isNovel','grp','mouse'};
[swrPtp_FN_lme_table] = get_lmetable([ripPtp(:,1); cleanRpRat(:,2)],cleanbvgrp,cleanMID,swrVarNames);
swrPtp_FN_lme = fitlme(swrPtp_FN_lme_table,'ripPtp ~ isNovel * grp + (1|mouse)');
swrPtp_FN_lmeF = plot_2wayLME(ripPtp(:,1)',ripPtp(:,2)',cleanbvgrp,lnlcols);
xticklabels({'F','N','F','N'}); ylim([0 1]); ylabel('P(SPW-R mod)')

if saveFlag
    fsave(swrRat_FN_lmeF,[sbase 'swr_lme_rate'])
    fsave(swrDur_FN_lmeF,[sbase 'swr_lme_duration'])
    fsave(swrLng_FN_lmeF,[sbase 'swr_lme_pLongRipple'])
    fsave(swrPtp_FN_lmeF,[sbase 'swr_lme_ripParticipation'])
end

%% SPWR Correlations to behavior
cleanLapRwd = uLapRwd;
cleanLapRwd(rmInds,:) = [];
cleanPsv = uPsv;
cleanPsv(rmInds,:) = [];

[psvXripPtp.mdl, psvXripPtp.fig] = plot_bhvXneur_corr(cleanPsv(:,1), ripPtp(:,1), cleanbvgrp, lnlcols);
xlabel('F LDI'); ylabel('F P(SPW-R mod'); ylim([0 1]); xlim([0 1])

[psvXripPtp.mdl, psvXripPtp.fig] = plot_bhvXneur_corr(cleanPsv(:,2), ripPtp(:,2), cleanbvgrp, lnlcols);
xlabel('N LDI'); ylabel('N P(SPW-R mod'); ylim([0 1]); xlim([-0.7 0.7])

[psvXripPtp.mdl, psvXripPtp.fig] = plot_bhvXneur_corr(cleanLapRwd(:,1), ripPtp(:,1), cleanbvgrp, lnlcols);
xlabel('F P(Rwd)'); ylabel('F P(SPW-R mod'); ylim([0 1]); xlim([0 1])

[psvXripPtp.mdl, psvXripPtp.fig] = plot_bhvXneur_corr(cleanLapRwd(:,2), ripPtp(:,2), cleanbvgrp, lnlcols);
xlabel('N P(Rwd)'); ylabel('N P(SPW-R mod'); ylim([0 1]); xlim([0 1])

%% Waterfall by sharp wave ripple peak
binedges = -wlen:histoBnsz:wlen;
nBins = length(binedges)-1;
ca1RipPk = find(binedges == 0,1);

% [swrBothPreSortPreFig,preMap,sortPre] = plot_unitWaterfall(rpMap(swrBothID,1:length(binedges)-1),binedges);
[swrBothPreSortPreFig,preMap,sortPre] = plot_unitWaterfall(rpMapZ(swrBothID,1:length(binedges)-1),binedges,0,1,0,0);
xticks([1,round(length(binedges)/2),nBins]); xticklabels(binedges([1,round(length(binedges)/2),length(binedges)]));
plot([ca1RipPk ca1RipPk],[0 sum(swrBothID)],'r--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Time to SWR peak (ms)')
[swrBothPreSortPreHisto,pkMapPre] = plot_unitPkHisto(preMap,binedges); xlim([-wlen wlen])
[~, ~, prePkBins] = pkChi2(pkMapPre,binedges);
xlabel('Time to SWR peak (ms)')

[swrBothPstSortPstFig,pstMap] = plot_unitWaterfall(rpMap(swrBothID,length(binedges):end),binedges);
plot([ca1RipPk ca1RipPk],[0 sum(swrBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Novel'); xlabel('Time to SWR peak (ms)')
[swrBothPstSortPstHisto,pkMapPst] = plot_unitPkHisto(pstMap,binedges); xlim([-wlen wlen])
[~, ~, pstPkBins] = pkChi2(pkMapPst,binedges);
xlabel('Time to SWR peak (ms)')

% [swrBothPstSortPreFig,pstMap] = plot_unitWaterfall(rpMap(swrBothID,length(binedges):end),binedges,sortPre,1,0);
[swrBothPstSortPreFig,pstMap] = plot_unitWaterfall(rpMapZ(swrBothID,length(binedges):end),binedges,sortPre,1,0,0);
xticks([1,round(length(binedges)/2),nBins]); xticklabels(binedges([1,round(length(binedges)/2),length(binedges)]));
plot([ca1RipPk ca1RipPk],[0 sum(swrBothID)],'r--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Time to SWR peak (ms)')

% Plot overlay pre to post
[~,ps.swr_PPPk_both,~,stats.swr_PPPk_both] = ttest(prePkBins,pstPkBins);
swrBothPkDistroHisto = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.14])
plot(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(pkMapPre)./sum(pkMapPre,'all'),'Color',fvncols(1,:))
plot(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(pkMapPst)./sum(pkMapPst,'all'),'Color',fvncols(2,:))
set(gca,'Position',[0.11 0.17 0.8 0.80]); xlim([-150 150])
xticks([binedges(1), binedges(round(nBins/2)+1), binedges(nBins+1)]); 
set(gca,'FontSize',12,'FontName','Arial')
text2bar(swrBothPkDistroHisto,'',ps.swr_PPPk_both,0.9);

swrCbar = plotColorbar(round([prctile(rpMapZ(swrBothID,length(binedges):end),1,'all'), prctile(rpMapZ(swrBothID,length(binedges):end),98,'all')],1),'parula');

if saveFlag
    fsave(swrBothPreSortPreFig,[sbase 'swr_both_pre_sortPre'])
    fsave(swrBothPstSortPstFig,[sbase 'swr_both_pst_sortPst'])
    fsave(swrBothPstSortPreFig,[sbase 'swr_both_pst_sortPre'])
    fsave(swrBothPreSortPreHisto,[sbase 'swr_both_pre_distro'])
    fsave(swrBothPstSortPstHisto,[sbase 'swr_both_pst_distro'])
    fsave(swrBothPkDistroHisto,[sbase 'swr_both_prepst_distro'])
    fsave(swrCbar,[sbase 'swr_Z_Cbar'])
end

%% Simple decoding analysis
% dcDat P(err) density FxF, NxN, NxF Mean (abs(err)) FxF, NxN, NxF

dcMat = [vertcat(dcDat.errLocsSub), vertcat(dcDat.errMeansSub)];
[~,ps.bdc_tt_FxF_lnnl,~,stats.bdc_tt_FxF_lnnl] = ttest2(dcMat(bvgrp(1).bvInd,4),dcMat(bvgrp(2).bvInd,4));
[~,ps.bdc_tt_NxN_lnnl,~,stats.bdc_tt_NxN_lnnl] = ttest2(dcMat(bvgrp(1).bvInd,5),dcMat(bvgrp(2).bvInd,5));
[~,ps.bdc_tt_NxF_lnnl,~,stats.bdc_tt_NxF_lnnl] = ttest2(dcMat(bvgrp(1).bvInd,6),dcMat(bvgrp(2).bvInd,6));

fxf_barF = plotBar2(dcMat(bvgrp(1).bvInd,4),dcMat(bvgrp(2).bvInd,4),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(fxf_barF, 'Mean Decode Error (m)', ps.bdc_tt_FxF_lnnl);
nxn_barF = plotBar2(dcMat(bvgrp(1).bvInd,5),dcMat(bvgrp(2).bvInd,5),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(nxn_barF, 'Mean Decode Error (m)', ps.bdc_tt_NxN_lnnl);
nxf_barF = plotBar2(dcMat(bvgrp(1).bvInd,6),dcMat(bvgrp(2).bvInd,6),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(nxf_barF, 'Mean Decode Error (m)', ps.bdc_tt_NxF_lnnl);

[~,ps.bdc_tt_FxF_lnnl_loc,~,stats.bdc_tt_FxF_lnnl_loc] = ttest2(dcMat(bvgrp(1).bvInd,1),dcMat(bvgrp(2).bvInd,1));
[~,ps.bdc_tt_NxN_lnnl_loc,~,stats.bdc_tt_NxN_lnnl_loc] = ttest2(dcMat(bvgrp(1).bvInd,2),dcMat(bvgrp(2).bvInd,2));
[~,ps.bdc_tt_NxF_lnnl_loc,~,stats.bdc_tt_NxF_lnnl_loc] = ttest2(dcMat(bvgrp(1).bvInd,3),dcMat(bvgrp(2).bvInd,3));

fxf_locF = plotBar2(dcMat(bvgrp(1).bvInd,1),dcMat(bvgrp(2).bvInd,1),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(fxf_locF, 'Error peak location (m)', ps.bdc_tt_FxF_lnnl_loc);
nxn_locF = plotBar2(dcMat(bvgrp(1).bvInd,2),dcMat(bvgrp(2).bvInd,2),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(nxn_locF, 'Error peak location (m)', ps.bdc_tt_NxN_lnnl_loc);
nxf_locF = plotBar2(dcMat(bvgrp(1).bvInd,3),dcMat(bvgrp(2).bvInd,3),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(nxf_locF, 'Error peak location (m)', ps.bdc_tt_NxF_lnnl_loc);

if saveFlag
    fsave(fxf_barF,[sbase 'decode_fxf_lnnl'],1,0);
    fsave(nxn_barF,[sbase 'decode_nxn_lnnl'],1,0);
    fsave(nxf_barF,[sbase 'decode_nxf_lnnl'],1,0);
    fsave(fxf_locF,[sbase 'decode_fxf_lnnl_loc'],1,0);
    fsave(nxn_locF,[sbase 'decode_nxn_lnnl_loc'],1,0);
    fsave(nxf_locF,[sbase 'decode_nxf_lnnl_loc'],1,0);
end

%% Correlating Decoding against LDI in each condition

preErr = [];
pstErr = [];
preLDI = [];
pstLDI = [];
mdlStats = [];
errXlap = [];

for i = 1:nMice
    cd(datT.fpath{i})

    epochfile = dir("*_dat.mat");
    load(epochfile.name)

    if mID(i) == 29 % Very crude editing to align trials
        if length(bvDat(i).pstLckDI) ~= length(dcDat(i).sub_nxn_lapAbsErr(sessLast.valTrials))
            bvDat(i).pstLckDI(end) = [];
        end
    end

    preLDI = [preLDI; bvDat(i).preLckDI];
    pstLDI = [pstLDI; bvDat(i).pstLckDI];
    preErr = [preErr; dcDat(i).sub_fxf_lapAbsErr(sessFrst.valTrials)];
    pstErr = [pstErr; dcDat(i).sub_nxn_lapAbsErr(sessLast.valTrials) dcDat(i).sub_nxf_lapAbsErr(sessLast.valTrials)];

    mdlStats(i).fxfRwdErr = mean(dcDat(i).sub_fxf_lapAbsErr(bvDat(i).preLapRwd),'omitnan');
    mdlStats(i).fxfNonErr = mean(dcDat(i).sub_fxf_lapAbsErr(~bvDat(i).preLapRwd),'omitnan');
    mdlStats(i).nxnRwdErr = mean(dcDat(i).sub_nxn_lapAbsErr(bvDat(i).pstLapRwd),'omitnan');
    mdlStats(i).nxnNonErr = mean(dcDat(i).sub_nxn_lapAbsErr(~bvDat(i).pstLapRwd),'omitnan');
    mdlStats(i).nxfRwdErr = mean(dcDat(i).sub_nxf_lapAbsErr(bvDat(i).pstLapRwd),'omitnan');
    mdlStats(i).nxfNonErr = mean(dcDat(i).sub_nxf_lapAbsErr(~bvDat(i).pstLapRwd),'omitnan');

    mdlFxF = get_linfit(bvDat(i).preLckDI, dcDat(i).sub_fxf_lapAbsErr(sessFrst.valTrials));
    mdlNxN = get_linfit(bvDat(i).pstLckDI, dcDat(i).sub_nxn_lapAbsErr(sessLast.valTrials));
    mdlNxF = get_linfit(bvDat(i).pstLckDI, dcDat(i).sub_nxf_lapAbsErr(sessLast.valTrials));
    [~, minindFxF] = min(bvDat(i).preLckDI);
    [~, minindNxN] = min(bvDat(i).pstLckDI);
    [~, maxindFxF] = max(bvDat(i).preLckDI);
    [~, maxindNxN] = max(bvDat(i).pstLckDI);

    mdlStats(i).fxfR = mdlFxF.r;
    mdlStats(i).nxnR = mdlNxN.r;
    mdlStats(i).nxfR = mdlNxF.r;
    mdlStats(i).fxfB = mdlFxF.b;
    mdlStats(i).nxnB = mdlNxN.b;
    mdlStats(i).nxfB = mdlNxF.b;

    try
        ldiXerrF = figure; hold on
        plot(bvDat(i).preLckDI, dcDat(i).sub_fxf_lapAbsErr(sessFrst.valTrials),'ko')
        plot(bvDat(i).pstLckDI, dcDat(i).sub_nxn_lapAbsErr(sessLast.valTrials),'ro')
        plot(bvDat(i).pstLckDI, dcDat(i).sub_nxf_lapAbsErr(sessLast.valTrials),'co')
        plot(bvDat(i).preLckDI([minindFxF maxindFxF]),mdlFxF.ypred([minindFxF maxindFxF]),'k');
        plot(bvDat(i).pstLckDI([minindNxN maxindNxN]),mdlNxN.ypred([minindNxN maxindNxN]),'r');
        plot(bvDat(i).pstLckDI([minindNxN maxindNxN]),mdlNxF.ypred([minindNxN maxindNxN]),'c');
        xlabel('Lap LDI'); ylabel('Lap Decoding Error (m)'); ylim([0 1]); xlim([-1 1]);
        set(gca,'FontName','Arial','FontSize',16)

        fsave(ldiXerrF, [rootFrst.name '_ldiXdecodeErr_velLo'], 1, 0);
    end
    close all
end
cd(groupSDir)

%% All trials Decoding against LDI Linear fit

mdlFxF = get_linfit(preLDI, preErr);
mdlNxN = get_linfit(pstLDI, pstErr(:,1));
mdlNxF = get_linfit(pstLDI, pstErr(:,2));

[~, minindFxF] = min(preLDI);
[~, minindNxN] = min(pstLDI);
[~, maxindFxF] = max(preLDI);
[~, maxindNxN] = max(pstLDI);

ldiXerrFxF = figure; hold on
plot(preLDI, preErr,'ko')
plot(pstLDI, pstErr(:,1),'ro')
plot(preLDI([minindFxF maxindFxF]),mdlFxF.ypred([minindFxF maxindFxF]),'k');
plot(pstLDI([minindNxN maxindNxN]),mdlNxN.ypred([minindNxN maxindNxN]),'r');
xlabel('Lap LDI'); ylim([0 1]); xlim([-1 1]);
set(gca,'FontName','Arial','FontSize',16)
text2bar(ldiXerrFxF, 'Lap Decoding Error (m)', mdlFxF.p, 0.7, 0.2, 'k');
text2bar(ldiXerrFxF, 'Lap Decoding Error (m)', mdlNxN.p, 0.4, 0.2, 'r');

ldiXerrNxF = figure; hold on
plot(pstLDI, pstErr(:,2),'co')
plot(pstLDI([minindNxN maxindNxN]),mdlNxF.ypred([minindNxN maxindNxN]),'c');
xlabel('Lap LDI'); ylim([0 1]); xlim([-1 1]);
set(gca,'FontName','Arial','FontSize',16)
text2bar(ldiXerrNxF, 'Lap Decoding Error (m)', mdlNxF.p, 0.6, 0.1, 'c');

% Compare Model R's Weak vs Strong
[~,ps.bdc_tt_fxfR,~,stats.bdc_tt_fxfR] = ttest2(vertcat(mdlStats(bvgrp(1).bvInd).fxfR),vertcat(mdlStats(bvgrp(2).bvInd).fxfR));
[~,ps.bdc_tt_nxnR,~,stats.bdc_tt_nxnR] = ttest2(vertcat(mdlStats(bvgrp(1).bvInd).nxnR),vertcat(mdlStats(bvgrp(2).bvInd).nxnR));
[~,ps.bdc_tt_nxfR,~,stats.bdc_tt_nxfR] = ttest2(vertcat(mdlStats(bvgrp(1).bvInd).nxfR),vertcat(mdlStats(bvgrp(2).bvInd).nxfR));
fxf_R_F = plotBar2(vertcat(mdlStats(bvgrp(1).bvInd).fxfR),vertcat(mdlStats(bvgrp(2).bvInd).fxfR),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(fxf_R_F, 'Pearson R', ps.bdc_tt_fxfR);
nxn_R_F = plotBar2(vertcat(mdlStats(bvgrp(1).bvInd).nxnR),vertcat(mdlStats(bvgrp(2).bvInd).nxnR),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(nxn_R_F, 'Pearson R', ps.bdc_tt_nxnR);
nxf_R_F = plotBar2(vertcat(mdlStats(bvgrp(1).bvInd).nxfR),vertcat(mdlStats(bvgrp(2).bvInd).nxfR),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([0 1]); text2bar(nxf_R_F, 'Pearson R', ps.bdc_tt_nxfR);

% Compare Model B's Weak vs Strong
[~,ps.bdc_tt_fxfB,~,stats.bdc_tt_fxfB] = ttest2(vertcat(mdlStats(bvgrp(1).bvInd).fxfB),vertcat(mdlStats(bvgrp(2).bvInd).fxfB));
[~,ps.bdc_tt_nxnB,~,stats.bdc_tt_nxnB] = ttest2(vertcat(mdlStats(bvgrp(1).bvInd).nxnB),vertcat(mdlStats(bvgrp(2).bvInd).nxnB));
[~,ps.bdc_tt_nxfB,~,stats.bdc_tt_nxfB] = ttest2(vertcat(mdlStats(bvgrp(1).bvInd).nxfB),vertcat(mdlStats(bvgrp(2).bvInd).nxfB));
fxf_B_F = plotBar2(vertcat(mdlStats(bvgrp(1).bvInd).fxfB),vertcat(mdlStats(bvgrp(2).bvInd).fxfB),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([-.5 .5]); text2bar(fxf_B_F, 'Lin. Model Slope', ps.bdc_tt_fxfB);
nxn_B_F = plotBar2(vertcat(mdlStats(bvgrp(1).bvInd).nxnB),vertcat(mdlStats(bvgrp(2).bvInd).nxnB),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([-.5 .5]); text2bar(nxn_B_F, 'Lin. Model Slope', ps.bdc_tt_nxnB);
nxf_B_F = plotBar2(vertcat(mdlStats(bvgrp(1).bvInd).nxfB),vertcat(mdlStats(bvgrp(2).bvInd).nxfB),lnlcols);
xticklabels({'Strong', 'Weak'}); ylim([-.5 .5]); text2bar(nxf_B_F, 'Lin. Model Slope', ps.bdc_tt_nxfB);

if saveFlag
    fsave(ldiXerrFxF, [sbase 'decode_fxf_ldiXerr'],1,0)
    fsave(ldiXerrNxF, [sbase 'decode_nxf_ldiXerr'],1,0)
    fsave(fxf_R_F, [sbase 'decode_fxf_linfit_r_lnnl'],1,0)
    fsave(nxn_R_F, [sbase 'decode_nxn_linfit_r_lnnl'],1,0)
    fsave(nxf_R_F, [sbase 'decode_nxf_linfit_r_lnnl'],1,0)
    fsave(fxf_B_F, [sbase 'decode_fxf_linfit_b_lnnl'],1,0)
    fsave(nxn_B_F, [sbase 'decode_nxn_linfit_b_lnnl'],1,0)
    fsave(nxf_B_F, [sbase 'decode_nxf_linfit_b_lnnl'],1,0)
end

%% Decoding error across laps

xsF = 1:50;
xsN = 51:100;

for i = 1:nMice
    lap_dcsub_fxf(i,xsF) = dcDat(i).sub_fxf_lapAbsErr(xsF+1); % Account for non-valid first lap per session
    % lap_dcsub_nxn(i,xsF) = dcDat(i).sub_nxn_lapAbsErr(xsF);
    lap_dcsub_nxf(i,xsF) = dcDat(i).sub_nxf_lapAbsErr(xsF);
end

mdl_dcsub_fxf_ln = get_linfit(xsF,mean(lap_dcsub_fxf(bvgrp(1).bvInd,xsF),'omitnan'));
mdl_dcsub_fxf_nl = get_linfit(xsF,mean(lap_dcsub_fxf(bvgrp(2).bvInd,xsF),'omitnan'));
% mdl_dcsub_nxn_ln = get_linfit(xsF,mean(lap_dcsub_nxn(bvgrp(1).bvInd,xsF),'omitnan'));
% mdl_dcsub_nxn_nl = get_linfit(xsF,mean(lap_dcsub_nxn(bvgrp(2).bvInd,xsF),'omitnan'));
mdl_dcsub_nxf_ln = get_linfit(xsF,mean(lap_dcsub_nxf(bvgrp(1).bvInd,xsF),'omitnan'));
mdl_dcsub_nxf_nl = get_linfit(xsF,mean(lap_dcsub_nxf(bvgrp(2).bvInd,xsF),'omitnan'));

[ciup_fxf_ln, cidn_fxf_ln] = get_CI(lap_dcsub_fxf(bvgrp(1).bvInd,xsF));
[ciup_fxf_nl, cidn_fxf_nl] = get_CI(lap_dcsub_fxf(bvgrp(2).bvInd,xsF));
% [ciup_nxn_ln, cidn_nxn_ln] = get_CI(lap_dcsub_nxn(bvgrp(1).bvInd,xsF));
% [ciup_nxn_nl, cidn_nxn_nl] = get_CI(lap_dcsub_nxn(bvgrp(2).bvInd,xsF));
[ciup_nxf_ln, cidn_nxf_ln] = get_CI(lap_dcsub_nxf(bvgrp(1).bvInd,xsF));
[ciup_nxf_nl, cidn_nxf_nl] = get_CI(lap_dcsub_nxf(bvgrp(2).bvInd,xsF));

lapErrLnF = figure; hold on
set(gcf,'Units','normalized','Position',[1.2 0.4 0.6315 0.1703])
plot_CIs(xsF, ciup_fxf_ln, cidn_fxf_ln, lnlcols(1,:)/2);
% plot_CIs(xsN, ciup_nxn_ln, cidn_nxn_ln, lnlcols(1,:));
plot_CIs(xsN, ciup_nxf_ln, cidn_nxf_ln, lnlcols(1,:)*1.5);
plot(xsF, mean(lap_dcsub_fxf(bvgrp(1).bvInd,xsF),'omitnan'), 'color', lnlcols(1,:)/2)
% plot(xsN, mean(lap_dcsub_nxn(bvgrp(1).bvInd,xsF),'omitnan'), 'color', lnlcols(1,:))
plot(xsN, mean(lap_dcsub_nxf(bvgrp(1).bvInd,xsF),'omitnan'), 'color', lnlcols(1,:)*1.5)
plot(xsF,mdl_dcsub_fxf_ln.ypred, 'color', lnlcols(1,:)/2, 'LineWidth',2)
% plot(xsN,mdl_dcsub_nxn_ln.ypred, 'color', lnlcols(1,:),   'LineWidth',2)
plot(xsN,mdl_dcsub_nxf_ln.ypred, 'color', lnlcols(1,:)*1.5,   'LineWidth',2)
xlabel("Lap #"); ylim([0 1])
text2bar(lapErrLnF,"", mdl_dcsub_fxf_ln.p, 0.8, 0.3, lnlcols(1,:)/2);
% text2bar(lapErrLnF,"", mdl_dcsub_nxn_ln.p, 0.3, 0.9, lnlcols(1,:));
text2bar(lapErrLnF,"Error (m)", mdl_dcsub_nxf_ln.p, 0.3, 0.1, lnlcols(1,:)*1.5);
set(gca,'FontSize',16,'FontName','Arial')

lapErrNlF = figure; hold on
set(gcf,'Units','normalized','Position',[1.2 0.4 0.6315 0.1703])
plot_CIs(xsF, ciup_fxf_nl, cidn_fxf_nl, lnlcols(2,:)/2);
% plot_CIs(xsN, ciup_nxn_nl, cidn_nxn_nl, lnlcols(2,:));
plot_CIs(xsN, ciup_nxf_nl, cidn_nxf_nl, lnlcols(2,:)*1.25);
plot(xsF, mean(lap_dcsub_fxf(bvgrp(2).bvInd,xsF),'omitnan'), 'color', lnlcols(2,:)/2)
% plot(xsN, mean(lap_dcsub_nxn(bvgrp(2).bvInd,xsF),'omitnan'), 'color', lnlcols(2,:))
plot(xsN, mean(lap_dcsub_nxf(bvgrp(2).bvInd,xsF),'omitnan'), 'color', lnlcols(2,:)*1.25)
plot(xsF,mdl_dcsub_fxf_nl.ypred, 'color', lnlcols(2,:)/2, 'LineWidth',2)
% plot(xsN,mdl_dcsub_nxn_nl.ypred, 'color', lnlcols(2,:),   'LineWidth',2)
plot(xsN,mdl_dcsub_nxf_nl.ypred, 'color', lnlcols(2,:)*1.25,   'LineWidth',2)
xlabel("Lap #"); ylim([0 1])
text2bar(lapErrNlF,"", mdl_dcsub_fxf_nl.p, 0.8, 0.3, lnlcols(2,:)/2);
% text2bar(lapErrNlF,"", mdl_dcsub_nxn_nl.p, 0.3, 0.7, lnlcols(2,:));
text2bar(lapErrNlF,"Error (m)", mdl_dcsub_nxf_nl.p, 0.3, 0.2, lnlcols(2,:)*1.25);
set(gca,'FontSize',16,'FontName','Arial')

% [~,ps.bdc_tt_nxn_ln_5v5,~,stats.bdc_tt_nxn_ln_5v5] = ttest2(mean(lap_dcsub_nxf(bvgrp(1).bvInd,1:5),2,'omitnan'),mean(lap_dcsub_nxf(bvgrp(1).bvInd,26:30),2,'omitnan'));
% plotBar2(mean(lap_dcsub_nxf(bvgrp(1).bvInd,1:5),2,'omitnan'),mean(lap_dcsub_nxf(bvgrp(1).bvInd,26:30),2,'omitnan'),lnlcols);

if saveFlag
    fsave(lapErrLnF,[sbase 'decode_errXlap_ln'],1,0);
    fsave(lapErrNlF,[sbase 'decode_errXlap_nl'],1,0);
end

%% Decoding error successful trials vs non-successful

dcErrVarNames = {'err','lapSuccess','grp','mouse'};
[dcErr_fxf_lme_table] = get_lmetable([vertcat(mdlStats.fxfRwdErr); vertcat(mdlStats.fxfNonErr)],bvgrp,mID,dcErrVarNames);
dcErr_fxf_lme = fitlme(dcErr_fxf_lme_table,'err ~ lapSuccess * grp + (1|mouse)');
dcErr_fxf_lmeF = plot_2wayLME(vertcat(mdlStats.fxfRwdErr)', vertcat(mdlStats.fxfNonErr)',bvgrp,lnlcols);
xticklabels({'Rwd','Miss','Rwd','Miss'}); ylabel('Error (m)'); ylim([0 1])

[dcErr_nxn_lme_table] = get_lmetable([vertcat(mdlStats.nxnRwdErr); vertcat(mdlStats.nxnNonErr)],bvgrp,mID,dcErrVarNames);
dcErr_nxn_lme = fitlme(dcErr_nxn_lme_table,'err ~ lapSuccess * grp + (1|mouse)');
dcErr_nxn_lmeF = plot_2wayLME(vertcat(mdlStats.nxnRwdErr)', vertcat(mdlStats.nxnNonErr)',bvgrp,lnlcols);
xticklabels({'Rwd','Miss','Rwd','Miss'}); ylabel('Error (m)'); ylim([0 1])

[dcErr_nxf_lme_table] = get_lmetable([vertcat(mdlStats.nxfRwdErr); vertcat(mdlStats.nxfNonErr)],bvgrp,mID,dcErrVarNames);
dcErr_nxf_lme = fitlme(dcErr_nxf_lme_table,'err ~ lapSuccess * grp + (1|mouse)');
dcErr_nxf_lmeF = plot_2wayLME(vertcat(mdlStats.nxfRwdErr)', vertcat(mdlStats.nxfNonErr)',bvgrp,lnlcols);
xticklabels({'Rwd','Miss','Rwd','Miss'}); ylabel('Error (m)'); ylim([0 1])

if saveFlag
    fsave(dcErr_fxf_lmeF,[sbase 'decode_fxf_rwd_lme'],1,0);
    fsave(dcErr_nxn_lmeF,[sbase 'decode_nxn_rwd_lme'],1,0);
    fsave(dcErr_nxf_lmeF,[sbase 'decode_nxf_rwd_lme'],1,0);
end

%% Compare Ripples across tr, ir, and rr cells

% SPWR Participation
[~,ps.swr_PPPre_trrr,~,stats.swr_PPPre_trrr] = ttest2(rpDat(trCells,2),rpDat(rrCells,2));
[~,ps.swr_PPPst_trrr,~,stats.swr_PPPst_trrr] = ttest2(rpDat(trCells,4),rpDat(rrCells,4));
[~,ps.swr_PPPrePst_tr,~,stats.swr_PPPrePst_tr] = ttest2(rpDat(trCells,2),rpDat(trCells,4));
[~,ps.swr_PPPrePst_rr,~,stats.swr_PPPrePst_rr] = ttest2(rpDat(rrCells,2),rpDat(rrCells,4));

swPartcpTRRRPreFig = plotBar2(rpDat(trCells,1),rpDat(rrCells,1));
xticklabels({'TR','RR'}); text2bar(swPartcpTRRRPreFig,'P(SPW-R Participation)',ps.swr_PPPre_trrr);
swPartcpTRRRPstFig = plotBar2(rpDat(trCells,3),rpDat(rrCells,3));
xticklabels({'TR','RR'}); text2bar(swPartcpTRRRPstFig,'P(SPW-R Participation)',ps.swr_PPPst_trrr);

bardat = histoBnsz*[mean(rpDat(trCells,2)), mean(rpDat(trCells,4)); mean(rpDat(rrCells,2)), mean(rpDat(rrCells,4))];
errdat = histoBnsz*[std(rpDat(trCells,2))./sqrt(sum(trCells)), std(rpDat(trCells,4))./sqrt(sum(trCells)); std(rpDat(rrCells,2))./sqrt(sum(rrCells)), std(rpDat(rrCells,4))./sqrt(sum(rrCells))];

swrModTRRRFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.12 0.27])
b1 = bar(bardat','FaceColor','flat');
b1(1).CData = [0 0 1];
b1(2).CData = [1 0 0];
errorbar([0.85 1.15 1.85 2.15],reshape(bardat,4,1),reshape(errdat,4,1),'k.')
xlim([0.5 2.5]); xticks(1:2); xticklabels({'Familiar', 'Novel'})
ylabel('Ripple Mod. Duration (ms)'); legend({'TR','RR'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

zmaxTRpre = max(abs(rpMapZ(trCells,1:length(binedges)-1)),[],2);
zmaxRRpre = max(abs(rpMapZ(rrCells,1:length(binedges)-1)),[],2);
zmaxTRpst = max(abs(rpMapZ(trCells,length(binedges):end)),[],2);
zmaxRRpst = max(abs(rpMapZ(rrCells,length(binedges):end)),[],2);

[~,ps.swr_ZPre_trrr,~,stats.swr_PPPre_trrr] = ttest2(zmaxTRpre,zmaxRRpre);
[~,ps.swr_ZPst_trrr,~,stats.swr_PPPst_trrr] = ttest2(zmaxTRpst,zmaxRRpst);
[~,ps.swr_ZPrePst_tr,~,stats.swr_PPPrePst_tr] = ttest2(zmaxTRpre, zmaxTRpst);
[~,ps.swr_ZPrePst_rr,~,stats.swr_PPPrePst_rr] = ttest2(zmaxRRpre, zmaxRRpst);

bardat = [mean(zmaxTRpre), mean(zmaxTRpst); mean(zmaxRRpre), mean(zmaxRRpst)];
errdat = [std(zmaxTRpre)./sqrt(sum(trCells)), std(zmaxTRpst)./sqrt(sum(trCells)); std(zmaxRRpre)./sqrt(sum(rrCells)), std(zmaxRRpst)./sqrt(sum(rrCells))];

% tdat = table([zeros(size(zmaxTRpre)); ones(size(zmaxRRpst))],[zmaxTRpre; zmaxRRpre],[zmaxTRpst; zmaxRRpst],...
% [rpDat(trCells,2); rpDat(rrCells,2)], [rpDat(trCells,4); rpDat(rrCells,4)],...
% 'VariableNames',{'unit_type','zPre','zPst','tPre','tPst'});
% tphase = table([1 2]','VariableNames',{'Phase'});
% rmZ = fitrm(tdat,'zPre-zPst~unit_type','WithinDesign',tphase);
% rmT = fitrm(tdat,'tPre-tPst~unit_type','WithinDesign',tphase);

swrZTRRRFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.12 0.27])
b1 = bar(bardat','FaceColor','flat');
b1(1).CData = [0 0 1];
b1(2).CData = [1 0 0];
errorbar([0.85 1.15 1.85 2.15],reshape(bardat,4,1),reshape(errdat,4,1),'k.')
xlim([0.5 2.5]); xticks(1:2); xticklabels({'Familiar', 'Novel'})
ylabel('Peak Z-score'); legend({'TR','RR'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

% % [~,ps.swr_postMod_trVrr,~,stats.swr_postMod_trVrr] = ttest(rpDat(:,1),rpDat(:,2));
% [~,ps.swr_prePrt_trVrr,~,stats.swr_prePrt_trVrr] = ttest2(rpDat(trCells,1),rpDat(rrCells,1));
% [~,ps.swr_pstPrt_trVrr,~,stats.swr_pstPrt_trVrr] = ttest2(rpDat(trCells,3),rpDat(rrCells,3));
% 
% swPartPreTRVRRFig = plotBar2(rpDat(trCells,1),rpDat(rrCells,1));
% ylabel('P(SPW-R Participation)'); xticklabels(["TR","RR"]); xlabel('Familiar')
% 
% swPartPstTRVRRFig = plotBar2(rpDat(trCells,3),rpDat(rrCells,3));
% ylabel('P(SPW-R Participation)'); xticklabels(["TR","RR"]); xlabel('Novel')
% 
% % swrTRRRCounts = [groupcounts(recID(trCells,1)) groupcounts(recID(rrCells,1))];
% % swrBothRatio = swrBothCounts ./ [groupcounts(recID(useCC,1)) groupcounts(recID(useCC,1))];
% % swrModCtFig = plot_barXmouse(swrBothRatio); ylim([0 1])
% % ylabel('P(Sig. Ripple Mod.)')

if saveFlag
    % saveas(swPartPreTRVRRFig,[sbase 'swr_trVrr_preParticip'],'png')
    % saveas(swPartPstTRVRRFig,[sbase 'swr_trVrr_pstParticip'],'png')
    fsave(swPartcpTRRRPreFig,[sbase 'swr_prePartbar_trrr'])
    fsave(swPartcpTRRRPstFig,[sbase 'swr_pstPartbar_trrr'])
    fsave(swrModTRRRFig,[sbase 'swr_ModDur_trrr'])
    fsave(swrZTRRRFig,[sbase 'swr_ModZ_trrr'])
end

%% Compare SPWRs against Behavior data
% swrBothCounts = [groupcounts(recID(swrFrstID,1)) groupcounts(recID(swrLastID,1))];
% uLckDI = [vertcat(bvDat.uPreLckDI), vertcat(bvDat.uPstLckDI)];

mdlLPreXspwr = get_linfit(uLDI(:,1),swrBothRatio(:,1));
mdlLPstXspwr = get_linfit(uLDI(:,2),swrBothRatio(:,2));

lckXspwrF = figure; hold on
plot(uLDI(:,1),swrBothRatio(:,1),'o','Color',fvncols(1,:))
plot(uLDI(:,2),swrBothRatio(:,2),'o','Color',fvncols(2,:))
plot(uLDI(:,1),mdlLPreXspwr.ypred,'Color',fvncols(1,:),'LineWidth',2)
plot(uLDI(:,2),mdlLPstXspwr.ypred,'Color',fvncols(2,:),'LineWidth',2)
xlim([-1 1]); xlabel('Lick DI');
ylim([0 1]); ylabel('P(SPWR-mod)');
legend({'Familiar','Novel'},'location','ne')
set(gca,'FontSize',16,'FontName','Arial')
xlims = xlim;
ylims = ylim;
text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims),  ['R = ' num2str(mdlLPreXspwr.r, 3)], 'Color', fvncols(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlLPreXspwr.p, 3)], 'Color', fvncols(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.25*diff(ylims), ['R = ' num2str(mdlLPstXspwr.r, 3)], 'Color', fvncols(2,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.3*diff(ylims),  ['p = ' num2str(mdlLPstXspwr.p, 3)], 'Color', fvncols(2,:), 'FontSize', 12)

% Compare delta of lckDI and delta of P(SPWR-Mod)
mdlDeltaLPreXspwr = get_linfit(uLDI(:,2)-uLDI(:,1),swrBothRatio(:,2)-swrBothRatio(:,1));
dlckXdspwrF = figure; hold on
plot(uLDI(:,2)-uLDI(:,1),swrBothRatio(:,2)-swrBothRatio(:,1),'ko')
plot(uLDI(:,2)-uLDI(:,1),mdlDeltaLPreXspwr.ypred,'Color','k','LineWidth',2)
xlabel('\Delta Lick DI');
ylabel('\Delta P(SPWR-mod)');
set(gca,'FontSize',16,'FontName','Arial')
xlims = xlim;ylims = ylim;
text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims),  ['R = ' num2str(mdlDeltaLPreXspwr.r, 3)], 'Color', fvncols(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlDeltaLPreXspwr.p, 3)], 'Color', fvncols(1,:), 'FontSize', 12)

if saveFlag
    fsave(lckXspwrF,[sbase 'swrXlck'])
    fsave(dlckXdspwrF,[sbase 'swrDXlckD'])
end

%% Velocity coding
% vlDat: 1&4 = sig.; 2&5 = Slope/B; 3&6 = Corr/R

% Summary of quantities
[velPieFig, vlprcts] = prepostPie(vlFrstID,vlLastID,useCC);
title("Velocity-modulated units");
vlBothCounts = [groupcounts(recID(vlFrstID,1)) groupcounts(recID(vlLastID,1))];
[~,ps.vl_ModCt_both,~,stats.vl_ModCt_both] = ttest(vlBothCounts(:,1),vlBothCounts(:,2));
vlBothCtFig = plot_barXmouse(vlBothCounts); ylim([0 70])
text2bar(vlBothCtFig,'# Sig. Velocity',ps.vl_ModCt_both);

[~,ps.vl_DtB_both,~,stats.vl_DtB_both] = ttest(vlDat(vlBothID,2), vlDat(vlBothID,5));
[~,ps.vl_DtR_both,~,stats.vl_DtR_both] = ttest(vlDat(vlBothID,3), vlDat(vlBothID,6));
[~,ps.vl_PPB_eith,~,stats.vl_PPB_eith] = ttest2(vlDat(vlFrstID,2),vlDat(vlLastID,5));
[~,ps.vl_PPR_eith,~,stats.vl_PPR_eith] = ttest2(vlDat(vlFrstID,3),vlDat(vlLastID,6));

% Units modulated in both phases
vlBFig = plotBar2(vlDat(vlBothID,2),vlDat(vlBothID,5)); ylim([-.5 1])
text2bar(vlBFig,'Velocity Slope (cm/s)',ps.vl_DtB_both);
vlRFig = plotBar2(vlDat(vlBothID,3),vlDat(vlBothID,6));
text2bar(vlRFig,'Velocity Corr. R^2',ps.vl_DtR_both);

% For units only modulated in one task phase
vlBEithFig = plotBar2(vlDat(vlFrstID & ~vlLastID,2),vlDat(vlLastID & ~vlFrstID,5)); ylim([-.5 1])
text2bar(vlBEithFig,'Velocity Slope (cm/s)',ps.vl_PPB_eith); xticklabels({'F-Only', 'N-Only'})
vlREithFig = plotBar2(vlDat(vlFrstID & ~vlLastID,3),vlDat(vlLastID & ~vlFrstID,6));
text2bar(vlREithFig,'Velocity Corr. R^2',ps.vl_PPR_eith); xticklabels({'F-Only', 'N-Only'})

% Delta histograms
vlDltB  = vlDat(vlBothID,5) - vlDat(vlBothID,2);
vlDltR  = vlDat(vlBothID,6) - vlDat(vlBothID,3);
binedges = -1:0.05:1;
vlDltRFig = plotDeltaHisto(vlDat(vlBothID,3),vlDat(vlBothID,6),binedges);
ylabel('Probability'); xlabel('\Delta R (Novel - Familiar)')
title("Significant units pre & post");
vlDltBFig = plotDeltaHisto(vlDat(vlBothID,2),vlDat(vlBothID,5),binedges);
xlabel('\Delta Slope (Novel - Familiar)')
title("Significant units pre & post");

if saveFlag
    fsave(velPieFig,[sbase 'vel_ModCt_pie'])
    fsave(vlBFig,[sbase 'vel_B'])
    fsave(vlRFig,[sbase 'vel_R'])
    fsave(vlBothCtFig,[sbase 'vel_ModCt_bar'])
    fsave(vlBEithFig,[sbase 'vel_B_eith'])
    fsave(vlREithFig,[sbase 'vel_R_eith'])
    % fsave(vlDltBFig,[sbase 'vel_deltaB'])
    % fsave(vlDltRFig,[sbase 'vel_deltaR'])
end

%% Theta Modulation
% thDat: 1&4 = sig.; 2&5 = MRL; 3&6 = Angle

[~,ps.th_PPM_both,~,stats.th_PPM_both] = ttest(thDat(thBothID,2),thDat(thBothID,5));
[~,ps.th_PPA_both,~,stats.th_PPA_both] = ttest(thDat(thBothID,3),thDat(thBothID,6));
[~,ps.th_PPM_eith,~,stats.th_PPM_eith] = ttest2(thDat(thFrstID & ~thLastID,2),thDat(thLastID & ~thFrstID,5));
[~,ps.th_PPA_eith,~,stats.th_PPA_eith] = ttest2(thDat(thFrstID & ~thLastID,3),thDat(thLastID & ~thFrstID,6));

% Summary of quantities
[thPieFig, thprcts] = prepostPie(thFrstID,thLastID,useCC);
title("Sig. Theta-modulated units");
thBothCounts = [groupcounts(recID(thFrstID,1)) groupcounts(recID(thLastID,1))];
[~,ps.th_ModCt_both,~,stats.th_ModCt_both] = ttest(thBothCounts(:,1),thBothCounts(:,2));
thBothCtFig = plot_barXmouse(thBothCounts); ylim([0 60]);
text2bar(thBothCtFig,'# Sig. Theta-Mod',ps.th_ModCt_both); 

% Units modulated in both phases
thMRLFig = plotBar2(thDat(thBothID,2),thDat(thBothID,5));
text2bar(thMRLFig,'Theta Mean Resultant Length',ps.th_PPM_both);
% 0/360 = trough, 180/540 = peak
thAngFig = plotBar2(thDat(thBothID,3),thDat(thBothID,6)); ylim([-180 180])
text2bar(thAngFig,'Theta Angle',ps.th_PPA_both);

% Delta histograms
% binedges = -0.1:0.01:0.1;
% thDltMFig = plotDeltaHisto(thDat(thBothID,2),thDat(thBothID,5),binedges);
% xlabel('\Delta Theta MRL (Novel - Familiar)')
% title("Significant units pre & post");
% binedges = rad2deg(-pi/2:pi/36:pi/2);
% thDltAFig = plotDeltaHisto(thDat(thBothID,3),thDat(thBothID,6),binedges);
% xlabel('\Delta Theta Angle (Novel - Familiar)')
% title("Significant units pre & post");

% For units only modulated in one task phase
thMRLEithFig = plotBar2(thDat(thFrstID & ~thLastID,2),thDat(thLastID & ~thFrstID,5));
text2bar(thMRLEithFig,'Theta Mean Resultant Length',ps.th_PPM_eith); xticklabels({'F-Only', 'N-Only'})
thAngEithFig = plotBar2(thDat(thFrstID & ~thLastID,3),thDat(thLastID & ~thFrstID,6)); ylim([-180 180])   % Compare frstHalf theta to last half theta
text2bar(thAngEithFig,'Theta Angle',ps.th_PPA_eith); xticklabels({'F-Only', 'N-Only'})

if saveFlag
    fsave(thPieFig,[sbase 'th_Mod_pie'])
    fsave(thMRLFig,[sbase 'th_MRL_bar'])
    fsave(thAngFig,[sbase 'th_Ang_bar'])
    fsave(thBothCtFig,[sbase 'th_ModCt_bar'])
    fsave(thMRLEithFig,[sbase 'th_MRL_eith'])
    fsave(thAngEithFig,[sbase 'th_Ang_eith'])
    % fsave(thDltMFig,[sbase 'th_MRL_delta'])
    % fsave(thDltAFig,[sbase 'th_Ang_delta'])
end

%% Waterfall by theta phase
binedges = rad2deg(0:pi/18:2*pi);
thPeak = find(binedges == 180,1);

[thBothPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(thMap(thBothID,1:length(binedges)-1),binedges);
plot([thPeak thPeak],[0 sum(thBothID)],'w--','LineWidth',2)
title('Familiar RZ, sort Familiar'); xlabel('Theta Phase')
thBothPreSortPreHisto = plot_unitPkHisto(tmpMap,binedges);
xlabel('Theta Phase')

[thBothPstSortPstFig,tmpMap] = plot_unitWaterfall(thMap(thBothID,length(binedges):end),binedges);
plot([thPeak thPeak],[0 sum(thBothID)],'w--','LineWidth',2)
title('Novel RZ, sort Novel'); xlabel('Theta Phase')
thBothPstSortPstHisto = plot_unitPkHisto(tmpMap,binedges);
xlabel('Theta Phase')

[thBothPstSortPreFig] = plot_unitWaterfall(thMap(thBothID,length(binedges):end),binedges,sortPre);
plot([thPeak thPeak],[0 sum(thBothID)],'w--','LineWidth',2)
title('Novel RZ, sort Familiar'); xlabel('Theta Phase')

if saveFlag
    saveas(thBothPreSortPreFig,[sbase 'th_both_pre_sortPre'],'png')
    saveas(thBothPstSortPstFig,[sbase 'th_both_pst_sortPst'],'png')
    saveas(thBothPstSortPreFig,[sbase 'th_both_pst_sortPre'],'png')
    saveas(thBothPreSortPreHisto,[sbase 'th_both_pre_distro'],'png')
    saveas(thBothPstSortPstHisto,[sbase 'th_both_pst_distro'],'png')
end

%% Theta phase figure
figure; hold on
cycleMax = 2*pi;
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.12])
% plot(rad2deg(0:pi/36:cycleMax),cos(pi:pi/36:5*pi),'k','LineWidth',2)    % Trough = 0/360
plot(rad2deg(0:pi/36:cycleMax),cos(0:pi/36:cycleMax),'k','LineWidth',2)     % Trough = 180
xlim([0 rad2deg(cycleMax)]); xticks([0 180 360 540 720]); xlabel('Theta Phase')
yticklabels(''); yticks([])
set(gca,'FontSize',12,'FontName','Arial')

%% Bar summary of SI, Vel, and Theta tuning

bardat = [siprcts; thprcts; vlprcts; swprcts];

figure; hold on; 
set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.39])
b = barh(bardat,'stacked');
b(1).FaceColor = [0.75 0.75 1];
b(2).FaceColor = fvncols(1,:);
b(3).FaceColor = fvncols(2,:);
b(4).FaceColor = [1 1 1];
yticks(1:4); yticklabels({'Spatial','Theta','Velocity','SPW-R'});
xlabel('Percentage Modulated'); ylim([0.5 4.5])
legend({'Both','Fam. only','Nov. only','Neither'})
box off
set(gca,'FontSize',16,'FontName','Arial')
%% Plot units by anatomical position and other data
proxCC = recID(:,5) < 0.5;
distCC = recID(:,5) > 0.5;
drslCC = recID(:,4) >= 0;
vtrlCC = recID(:,4) < 0;

[~,ps.fr_runn_proxdist,~,stats.fr_runn_proxdist] = ttest2(frDat(proxCC,2),frDat(distCC,2));
[~,ps.fr_stnd_proxdist,~,stats.fr_stnd_proxdist] = ttest2(frDat(proxCC,1),frDat(distCC,1));
[~,ps.fr_runn_drslvtrl,~,stats.fr_runn_drslvtrl] = ttest2(frDat(drslCC,2),frDat(vtrlCC,2));
[~,ps.fr_stnd_drslvtrl,~,stats.fr_stnd_drslvtrl] = ttest2(frDat(drslCC,1),frDat(vtrlCC,1));
[~,ps.si_proxdist,~,stats.si_proxdist]           = ttest2(lcDat(proxCC,2),lcDat(distCC,2));
[~,ps.si_drslvtrl,~,stats.si_drslvtrl]           = ttest2(lcDat(drslCC,2),lcDat(vtrlCC,2));
[~,ps.th_ang_proxdist,~,stats.th_ang_proxdist]   = ttest2(thDat(proxCC,3),thDat(distCC,3));
[~,ps.th_ang_drslvtrl,~,stats.th_ang_drslvtrl]   = ttest2(thDat(drslCC,3),thDat(vtrlCC,3));
[~,ps.swr_mod_proxdist,~,stats.swr_mod_proxdist] = ttest2(rpDat(proxCC,2),rpDat(distCC,2));
[~,ps.swr_mod_drslvtrl,~,stats.swr_mod_drslvtrl] = ttest2(rpDat(drslCC,4),rpDat(vtrlCC,4));

% Compare theta pre shift for significant units pre-shift only
thBothID = useCC & thDat(:,1) <= 0.05; % & thDat(:,4) <= 0.05;
pdCMap = make_custom_cmap([0.063 0.322 0.255], [0.353 0.612 0.223],sum(thBothID));
dvCMap = make_custom_cmap([0.353 0.835 0.772], [0.063 0.482 0.416],sum(thBothID));
datXanatThAngFig = plot_group_datXlyr(thDat(:,3)+180,thBothID,recID(:,4:5),hsv(sum(thBothID)));
xlim([0 0.8]); ylim([-200 250]); % datXanatThAngFig.Children(1).Label.String = 'Theta phase \circ'; clim([0 360])
% thCbar = plotColorbar([-180 180],'hsv');
% thAngProxDistFig = plot_proxVdist(thDat(:,3)+180,thBothID & proxCC,thBothID & distCC);
% text2bar(thAngProxDistFig,'Theta phase \circ',ps.th_ang_proxdist); ylim([0 360]);
% thAngDrslVtrlFig = plot_drslVvtrl(thDat(:,3)+180,thBothID & drslCC,thBothID & vtrlCC);
% xlabel('Theta phase \circ'); xlim([0 360]); text2bar(thAngDrslVtrlFig,'',ps.th_ang_drslvtrl);
[thAngPDCorrFig,mdlparams] = plot_anatCorr(thDat(thBothID,3)+180,recID(thBothID,5),1,pdCMap);
ylim([0 360]); text2corr(thAngPDCorrFig,'Theta phase \circ',mdlparams); 
[thAngDVCorrFig,mdlparams] = plot_anatCorr(thDat(thBothID,3)+180,recID(thBothID,4),2,dvCMap);
xlabel('Theta phase \circ'); xlim([0 360]); text2corr(thAngDVCorrFig,'',mdlparams);

%Compare run FR pre shift
pdCMap = make_custom_cmap([0.063 0.322 0.255], [0.353 0.612 0.223],sum(useCC));
dvCMap = make_custom_cmap([0.353 0.835 0.772], [0.063 0.482 0.416],sum(useCC));
datXanatFRFig = plot_group_datXlyr(frDat(:,2),useCC,recID(:,4:5),hot(sum(useCC)));
xlim([0 0.8]); ylim([-200 250]); % datXanatFRFig.Children(1).Label.String = 'Firing Rate (Hz)';
% frCbar = plotColorbar([0 round(max(frDat(useCC,2)))],'hot');
% frProxDistFig = plot_proxVdist(frDat(:,2),useCC & proxCC,useCC & distCC);
% text2bar(frProxDistFig,'Firing Rate (Hz)',ps.fr_runn_proxdist);
% frDrslVtrlFig = plot_drslVvtrl(frDat(:,2),useCC & drslCC,useCC & vtrlCC);
% xlabel('Firing Rate (Hz)'); text2bar(frDrslVtrlFig,'',ps.fr_runn_drslvtrl);
[frPVCorrFig,mdlparams] = plot_anatCorr(frDat(useCC,2),recID(useCC,5),1,pdCMap);
text2corr(frPVCorrFig,'Running FR (Hz)',mdlparams);
[frDVCorrFig,mdlparams] = plot_anatCorr(frDat(useCC,2),recID(useCC,4),2,dvCMap);
xlabel('Running FR (Hz)'); text2corr(frDVCorrFig,'',mdlparams);

% Compare SI pre shift for significant units pre-shift only
pdCMap = make_custom_cmap([0.063 0.322 0.255], [0.353 0.612 0.223],sum(siFrstID));
dvCMap = make_custom_cmap([0.353 0.835 0.772], [0.063 0.482 0.416],sum(siFrstID));
datXanatSIFig = plot_group_datXlyr(lcDat(:,2),siFrstID,recID(:,4:5),hot(sum(siFrstID)));
clim([0 prctile(lcDat(siFrstID,2),99)])
xlim([0 0.8]); ylim([-200 250]); % datXanatSIFig.Children(1).Label.String = 'Spatial Info. (bits/spike)';
% siCbar = plotColorbar([0 round(prctile(lcDat(siBothID,2),99))],'hot');
% siProxDistFig = plot_proxVdist(lcDat(:,2),siBothID & proxCC,siBothID & distCC); ylim([0 prctile(lcDat(siBothID,2),99)]);
% text2bar(siProxDistFig,'Spatial Info. (bits/spike)',ps.si_proxdist);
% siDrslVtrlFig = plot_drslVvtrl(lcDat(:,2),siBothID & drslCC,siBothID & vtrlCC);
% xlabel('Spatial Info. (bits/spike)'); text2bar(siDrslVtrlFig,'',ps.si_drslvtrl); xlim([0 prctile(lcDat(siBothID,2),99)]);
[siPDCorrFig,mdlparams] = plot_anatCorr(lcDat(siFrstID,2),recID(siFrstID,5),1,pdCMap);
ylim([0 4]); text2corr(siPDCorrFig,'Spatial Info. (bits/spike)',mdlparams);
[siDVCorrFig,mdlparams] = plot_anatCorr(lcDat(siFrstID,2),recID(siFrstID,4),2,dvCMap);
xlim([0 4]); xlabel('Spatial Info. (bits/spike)'); text2corr(siDVCorrFig,'',mdlparams);
%%
% Compare BI pre shift for significant units pre-shift only
pdCMap = make_custom_cmap([0.063 0.322 0.255], [0.353 0.612 0.223],sum(bstFrstID));
dvCMap = make_custom_cmap([0.353 0.835 0.772], [0.063 0.482 0.416],sum(bstFrstID));
datXanatBIFig = plot_group_datXlyr(bsDat(:,1),bstFrstID,recID(:,4:5),hot(sum(bstFrstID)));
clim([0 prctile(bsDat(bstFrstID,2),99)])
xlim([0 0.8]); ylim([-200 250]); % datXanatSIFig.Children(1).Label.String = 'Spatial Info. (bits/spike)';
% siCbar = plotColorbar([0 round(prctile(lcDat(siBothID,2),99))],'hot');
% siProxDistFig = plot_proxVdist(lcDat(:,2),siBothID & proxCC,siBothID & distCC); ylim([0 prctile(lcDat(siBothID,2),99)]);
% text2bar(siProxDistFig,'Spatial Info. (bits/spike)',ps.si_proxdist);
% siDrslVtrlFig = plot_drslVvtrl(lcDat(:,2),siBothID & drslCC,siBothID & vtrlCC);
% xlabel('Spatial Info. (bits/spike)'); text2bar(siDrslVtrlFig,'',ps.si_drslvtrl); xlim([0 prctile(lcDat(siBothID,2),99)]);
[biPDCorrFig,mdlparams] = plot_anatCorr(bsDat(bstFrstID,2),recID(bstFrstID,5),1,pdCMap);
ylim([0 4]); text2corr(biPDCorrFig,'Spatial Info. (bits/spike)',mdlparams);
[biDVCorrFig,mdlparams] = plot_anatCorr(bsDat(bstFrstID,2),recID(bstFrstID,4),2,dvCMap);
xlim([0 4]); xlabel('Spatial Info. (bits/spike)'); text2corr(biDVCorrFig,'',mdlparams);
%%
% % Compare Ripple Mod Dur pre shift for significant units pre-shift only
% swrBothID = useCC & rpDat(:,2) > 0; % & rpDat(:,5) <= 0.05;
% datXanatSWRFig = plot_group_datXlyr(rpDat(:,2)*histoBnsz,swrBothID,recID(:,4:5),turbo(sum(swrBothID)));
% clim([0 prctile(rpDat(swrBothID,2)*histoBnsz,95)])
% xlim([0 0.8]); ylim([-200 250]); % datXanatSIFig.Children(1).Label.String = 'Spatial Info. (bits/spike)';
% swrCbar = plotColorbar([0 prctile(rpDat(swrBothID,2)*histoBnsz,95)],'turbo');
% swrProxDistFig = plot_proxVdist(rpDat(:,2)*histoBnsz,swrBothID & proxCC,swrBothID & distCC); ylim([0 prctile(rpDat(swrBothID,2),95)]);
% text2bar(swrProxDistFig,'SWR Mod. Dur. (ms)',ps.swr_mod_proxdist);
% swrDrslVtrlFig = plot_drslVvtrl(rpDat(:,2)*histoBnsz,swrBothID & drslCC,swrBothID & vtrlCC);
% xlabel('SWR Mod. Dur. (ms)'); text2bar(swrDrslVtrlFig,'',ps.swr_mod_drslvtrl); xlim([0 prctile(rpDat(swrBothID,2),95)]);
% swrPDCorrFig = plot_anatCorr(rpDat(swrBothID,2)*histoBnsz,recID(swrBothID,5),1,winter(sum(swrBothID)));
% ylabel('SWR Mod. Dur. (ms)')
% swrDVCorrFig = plot_anatCorr(rpDat(swrBothID,2)*histoBnsz,recID(swrBothID,4),2,parula(sum(swrBothID)));
% xlabel('SWR Mod. Dur. (ms)')

typeCells = zeros(size(recID,1),1); typeCells(rrCells) = 1; typeCells(trCells) = -1;
datXanatRFFig = plot_group_datXlyr(typeCells,siBothID,recID(:,4:5),winter(sum(siBothID)));
xlim([0 0.8]); ylim([-200 250]); colormap(redbluecmap(3)); % datXanatRFFig.Children(1).Ticks = [0 1];
celltypeCbar = plotColorbar({'TR','RR'},redbluecmap(3));
% rfProxDistFig = plot_proxVdist(typeCells,siBothID & proxCC,siBothID & distCC);
% ylabel('Spatial Info. (bits/spike)')
% rfDrslVtrlFig = plot_drslVvtrl(typeCells,siBothID & drslCC,siBothID & vtrlCC);
% xlabel('Spatial Info. (bits/spike)')
% rfPDCorrFig = plot_anatCorr(typeCells(siBothID),recID(siBothID,5),1,winter(sum(siBothID)));
% ylabel('Spatial Info. (bits/spike)')
% rfDVCorrFig = plot_anatCorr(typeCells(siBothID),recID(siBothID,4),2,parula(sum(siBothID)));
% xlabel('Spatial Info. (bits/spike)')

if saveFlag
    fsave(datXanatThAngFig,[sbase 'anat_thAng_PDDV_pre'])
    % fsave(thAngProxDistFig,[sbase 'anat_thAng_PDbar_pre'])
    % fsave(thAngDrslVtrlFig,[sbase 'anat_thAng_DVbar_pre'])
    fsave(thAngPDCorrFig,[sbase 'anat_thAng_PDcor_pre'])
    fsave(thAngDVCorrFig,[sbase 'anat_thAng_DVcor_pre'])
    fsave(datXanatFRFig,[sbase 'anat_frRun_PDDV_pre'])
    % fsave(frProxDistFig,[sbase 'anat_frRun_PDbar_pre'])
    % fsave(frDrslVtrlFig,[sbase 'anat_frRun_DVbar_pre'])
    fsave(frPVCorrFig,[sbase 'anat_frRun_PDcor_pre'])
    fsave(frDVCorrFig,[sbase 'anat_frRun_DVcor_pre'])
    fsave(datXanatSIFig,[sbase 'anat_si_PDDV_pre'])
    % fsave(siProxDistFig,[sbase 'anat_si_PDbar_pre'])
    % fsave(siDrslVtrlFig,[sbase 'anat_si_DVbar_pre'])
    fsave(siPDCorrFig,[sbase 'anat_si_PDcor_pre'])
    fsave(siDVCorrFig,[sbase 'anat_si_DVcor_pre'])
    fsave(datXanatRFFig,[sbase 'anat_refFrame_PDDV'],'png')
    % fsave(frCbar,[sbase 'anat_frRun_cbar'])
    % fsave(thCbar,[sbase 'anat_thA_cbar'])
    % fsave(siCbar,[sbase 'anat_si_cbar'])
    fsave(celltypeCbar,[sbase 'anat_refFrame_cbar'])
    fsave(datXanatBIFig,[sbase 'anat_bi_PDDV_pre'])
    fsave(biPDCorrFig,[sbase 'anat_bi_PDcor_pre'])
    fsave(biDVCorrFig,[sbase 'anat_bi_DVcor_pre'])
end

%% Scatter TR and RR by depth/distance normalized to all spatially mod cell

dvbnsz = 10; dvBinEdges = min(recID(siBothID,4))-dvbnsz/2:dvbnsz:max(recID(siBothID,4))+0.5*dvbnsz; dvBinCtrs = dvBinEdges(2:end)-0.5*dvbnsz;
pdbnsz = 0.05; pdBinEdges = 0:pdbnsz:round(max(recID(siBothID,5)),1)+pdbnsz/2; pdBinCtrs = pdBinEdges(2:end)-0.5*pdbnsz;
trDat = recID(trCells,4:5);
rrDat = recID(rrCells,4:5);
irDat = recID(irCells,4:5);
allDat = recID(siBothID,4:5);

dvCountTR = histcounts(trDat(:,1),dvBinEdges);
dvCountRR = histcounts(rrDat(:,1),dvBinEdges);
dvCountIR = histcounts(irDat(:,1),dvBinEdges);
dvCountAll = histcounts(allDat(:,1),dvBinEdges);
pdCountTR = histcounts(trDat(:,2),pdBinEdges);
pdCountRR = histcounts(rrDat(:,2),pdBinEdges);
pdCountIR = histcounts(irDat(:,2),pdBinEdges);
pdCountAll = histcounts(allDat(:,2),pdBinEdges);

% Models
dvRRmdl = get_linfit(dvBinCtrs,dvCountRR);
dvTRmdl = get_linfit(dvBinCtrs,dvCountTR);
pdRRmdl = get_linfit(pdBinCtrs,pdCountRR);
pdTRmdl = get_linfit(pdBinCtrs,pdCountTR);

rfPDCorrFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
scatter(pdBinCtrs,pdCountRR./sum(pdCountAll,'all'),'r','MarkerFaceColor','r')
plot(pdBinCtrs,pdRRmdl.ypred./sum(pdCountAll,'all'),'r')
scatter(pdBinCtrs,pdCountTR./sum(pdCountAll,'all'),'b','MarkerFaceColor','b')
plot(pdBinCtrs,pdTRmdl.ypred./sum(pdCountAll,'all'),'b')
xlabel('% Distance through subiculum'); ylim([0 .2])
text2corr(rfDVCorrFig,'% of spatial cells',pdTRmdl,0.4);
text2corr(rfDVCorrFig,'% of spatial cells',pdRRmdl,0.9);
set(gca,'FontSize',12,'FontName','Arial')

rfDVCorrFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
scatter(dvCountRR./sum(dvCountAll,'all'),dvBinCtrs,'r','MarkerFaceColor','r')
plot(dvRRmdl.ypred./sum(dvCountAll,'all'),dvBinCtrs,'r')
scatter(dvCountTR./sum(dvCountAll,'all'),dvBinCtrs,'b','MarkerFaceColor','b')
plot(dvTRmdl.ypred./sum(dvCountAll,'all'),dvBinCtrs,'b')
xlabel('% of spatial cells')
set(gca,'FontSize',12,'FontName','Arial')
xlim([0 .1]); xlabel('Spatial Info. (bits/spike)'); 
text2corr(rfDVCorrFig,'Distance to layer center (um)',dvRRmdl,0.4);
text2corr(rfDVCorrFig,'Distance to layer center (um)',dvTRmdl,0.9);

if saveFlag
    fsave(rfPDCorrFig,[sbase 'anat_refFrame_PDcor'])
    fsave(rfDVCorrFig,[sbase 'anat_refFrame_DVcor'])
end

%%
% Count units by type and reference frame for each mouse and shank

% mID; shank; pos; good units; tr; rr; ir;
rfMat = [];
for i = 1:nMice
    units = recID(:,1) == mID(i);
    shs = unique(recID(units,5));
    for j = 1:length(shs)
        if shs(j) < 0
            continue
        end
        shUnits = recID(:,5) == shs(j) & recID(:,1) == mID(i);
        nTot = sum(shUnits & useCC);
        nTR = sum(shUnits & trCells);
        nRR = sum(shUnits & rrCells);
        nIR = sum(shUnits & irCells);
        rfMat = [rfMat; mID(i), j, shs(j), nTot, nTR, nRR, nIR];
    end
end

pTR = rfMat(:,5) ./ rfMat(:,4);
pRR = rfMat(:,6) ./ rfMat(:,4);
pIR = rfMat(:,7) ./ rfMat(:,4);
xDat = [rfMat(:,3); rfMat(:,3)];
yDat = [pTR; pRR];
cDat = [-1*ones(size(pTR)); ones(size(pRR))];

pdTRmdl = get_linfit(rfMat(:,3),pTR);
pdRRmdl = get_linfit(rfMat(:,3),pRR);

rfPDCorrFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
scatter(xDat,yDat,[],cDat,'filled')
colormap(redbluecmap(3))
plot(rfMat(~isnan(pTR),3),pdTRmdl.ypred,'b')
plot(rfMat(~isnan(pRR),3),pdRRmdl.ypred,'r')
xlabel('% Distance through subiculum'); ylim([0 1])
text2corr(rfPDCorrFig,'',pdTRmdl,0.4);
text2corr(rfPDCorrFig,'P(Unit type by shank)',pdRRmdl,0.9);
set(gca,'FontSize',16,'FontName','Arial')

if saveFlag
    fsave(rfPDCorrFig,[sbase 'anat_refFrame_PDcor'])
end

%% Save stats

if saveFlag
    save([sbase, 'stats'],'ps','stats')
end

%% Functions

function [fhandle,prcts] = prepostPie(preSig,postSig,useCCs)

bothSig = preSig & postSig;
notSig = not(preSig | postSig) & useCCs;

nBoth = sum(bothSig);
nPre = sum(preSig);
nPost = sum(postSig);
nNot = sum(notSig);

prcts = [nBoth, nPre - nBoth, nPost - nBoth, nNot] ./ sum(useCCs); 

cMap = [0.25 0.15 1; 0.5 0.5 1; 0.75 0.75 1; 0.25 0.25 0.25];

fhandle = figure;
p = piechart([nBoth, nPre - nBoth, nPost - nBoth, nNot],["Both","Familiar-only","Novel-only","Neither"]);
p.LabelStyle = 'namedata';
colororder(cMap)
end

function [fhandle] = plot_bhvrTraceCI(dat1,dat2,vcolors)

[ciup1, cidn1] = get_CI(dat1);
[ciup2, cidn2] = get_CI(dat2);

bnpos = linspace(-92.5,92.5,size(dat1,2));

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.5 0.45 0.14])

plot_CIs(bnpos,ciup1,cidn1,vcolors(1,:))
plot(bnpos,mean(dat1),'Color',vcolors(1,:),'LineWidth',2)
plot_CIs(bnpos,ciup2,cidn2,vcolors(2,:))
plot(bnpos,mean(dat2),'Color',vcolors(2,:),'LineWidth',2)
xlabel('Distance to RZ (cm)')
set(gca,'FontSize',16,'FontName','Arial')

end

function [fhandle] = plot_3bhvrTraceCI(dat,grpInds,vcolors)

nGrps = length(grpInds);

for i = 1:nGrps
    [cis(i).up, cis(i).dn] = get_CI(dat(grpInds(i).inds,:));
end

bnpos = linspace(-92.5,92.5,size(dat,2));

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.5 0.45 0.14])

for i = 1:nGrps
    plot_CIs(bnpos,cis(i).up,cis(i).dn,vcolors(i,:))
    plot(bnpos, mean(dat(grpInds(i).inds,:)),'Color',vcolors(i,:),'LineWidth',2)
end
xlabel('Distance to RZ (cm)')
set(gca,'FontSize',16,'FontName','Arial')

end

function [fhandle] = plot_multiBar(dat,grpInds,vcolors)

nGrps = length(grpInds);
grpID = zeros(size(dat));

for i = 1:nGrps
    grpMean(i) = mean(dat(grpInds(i).inds,1));
    grpsem(i) = std(dat(grpInds(i).inds,1)) ./ sqrt(length(grpInds(i).inds));
    grpID(grpInds(i).inds) = i;
end

fhandle = figure; hold on 
set(gcf,'units','normalized','position',[0.4 0.5 0.22 0.23])
b = bar(1:nGrps,grpMean,'FaceColor','flat');
b.CData = vcolors;
errorbar(1:nGrps,grpMean,grpsem,'k.')
plot(grpID,dat,'o','Color',[.7 .7 .7])
xlim([0.5 nGrps + 0.5])
set(gca,'FontSize',16,'FontName','Arial')

end

function [fhandle] = plotMiniBar(dat1,dat2,vColors)

arguments
    dat1
    dat2
    % vColors = [0.5 0.5 1; 0.75 0.75 1];
    vColors = [.35 .35 .35; 1 .25 .25];
end

nUnits = [size(dat1,1) size(dat2,1)];
bardat = [mean(dat1); mean(dat2)];
semdat = [std(dat1)/sqrt(nUnits(1)); std(dat2)/sqrt(nUnits(2))];
xrands1 = (rand(nUnits(1),1)-0.5)*0.2;
xrands2 = (rand(nUnits(2),1)-0.5)*0.2;

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.5 0.25 0.13])
% set(gcf,'units','normalized','position',[0.4 0.35 0.1 0.2])
% b = bar([1.15 2.15],bardat,0.3,'FaceColor','flat','BarWidth',0.5);
% b.CData = vColors;
plot(xrands1+1,dat1,'.','Color',vColors(1,:),'MarkerSize',10)
plot(xrands2+2,dat2,'.','Color',vColors(2,:),'MarkerSize',10)
errorbar([1.15 2.15],bardat,semdat,'k.')
xlim([0.5 2.5]); ylim([0 max([prctile(dat1,98); prctile(dat2,98)],[],'all')])
xticks(1:2); xticklabels({'F', 'N'})
set(gca,'FontSize',16,'FontName','Arial')
box off
end

function [fhandle] = plotPVCorrComp(dat1, dat2, nMice, pval)
% Plot on-diagonal vs off diagonal average PVC

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.2])
plot([dat1' dat2']','-o','Color',[.5 .5 .5])
errorbar([1 2],mean([dat1' dat2'],1,'omitnan'),std([dat1' dat2'],1,'omitnan')./sqrt(nMice),'k.','LineWidth',2,'CapSize',20)
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Diagonal','Off-Diag'})
ylim([-0.25 1]); text2bar(fhandle,'Mean of PV Corr.',pval);
set(gca,'FontSize',16,'FontName','Arial')

end

function [fhandle] = plotDistroHisto(distro1,distro2,binpos,rzPos)
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

fhandle = figure; hold on
plot(binpos*100,distro1./sum(distro1),'Color',vColors2(1,:),'LineWidth',2);
plot(binpos*100,distro2./sum(distro2),'Color',vColors2(2,:),'LineWidth',2);
if length(rzPos) == 2
    plot([rzPos(1) rzPos(1)]*100,[0 0.15],'--','Color',vColors2(1,:))
    plot([rzPos(2) rzPos(2)]*100,[0 0.15],'--','Color',vColors2(2,:))
else
    plot([rzPos rzPos]*100,[0 0.15],'k--')
end
ylim([0 max([distro1./sum(distro1); distro2./sum(distro2)+0.02],[],'all')])
xlim([binpos(1)*100-1 binpos(end)*100+1])
ylabel('P(Field Peak)')
legend({'Familiar','Novel'},'Location','northeast')
set(gca,'FontSize',12,'FontName','Arial')
end

function [fhandle] = plot_2d_bhvr(dat,learn,nonlearn,cols)
fhandle = figure; hold on; axis square
set(gcf,'Units','normalized','Position',[0.1 0.4 0.3333 0.1786])
plot(dat(learn,1),dat(learn,2),'o','Color',cols(1,:))
plot(dat(nonlearn,1),dat(nonlearn,2),'o','Color',cols(2,:))
set(gca,'FontSize',16,'FontName','Arial')
end

function [lmetbl] = get_lmetable(dat,bvgrp,mID,varnames)
% Assumes dat = [condition1; condition2] where length(condition1) = bvgrp(1).n

splitvar = [zeros(length(dat)/2,1); ones(length(dat)/2,1)];
expgrp = zeros(length(dat)/2,1);
expgrp(bvgrp(1).bvInd) = 1;
expgrp = repmat(expgrp,[2,1]);
lmetbl = table(dat,splitvar,expgrp,[mID; mID],'VariableNames',varnames);
end

function [fhandle] = plot_2wayLME(dat1,dat2,bvgrp,cols)
fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.2 0.2])
fhandle = fixRatio(fhandle);
for i = 1:2
    nMice = bvgrp(i).n;
    plot([0.85*ones(nMice,1) 1.15*ones(nMice,1)]'+(i-1), [dat1(bvgrp(i).bvInd)' dat2(bvgrp(i).bvInd)']','-o','Color',cols(i,:))
    errorbar([0.85 1.15]+(i-1),mean([dat1(bvgrp(i).bvInd)' dat2(bvgrp(i).bvInd)'],1,'omitnan'),std([dat1(bvgrp(i).bvInd)' dat2(bvgrp(i).bvInd)'],1,'omitnan')./sqrt(nMice),'k.','LineWidth',2,'CapSize',20)
end
xlim([0.5 2.5]); xticks([0.85 1.15 1.85 2.15]);
ylim([-0.57 1]);
set(gca,'FontSize',16,'FontName','Arial')
end

function [mdl, fhandle] = plot_bhvXneur_corr(xdat,ydat, cleanbvgrp, cols)
mdl = get_linfit(xdat,ydat);

fhandle = figure; hold on
plot(xdat,mdl.ypred,'k','LineWidth',1)
plot(xdat(cleanbvgrp(1).bvInd),ydat(cleanbvgrp(1).bvInd),'o','Color',cols(1,:));
plot(xdat(cleanbvgrp(2).bvInd),ydat(cleanbvgrp(2).bvInd),'o','Color',cols(2,:));
text2bar(fhandle, '', mdl.p);
set(gcf,'units','normalized','position',[0.4 0.35 0.4 0.14])
set(gca,'FontSize',16,'FontName','Arial')
end