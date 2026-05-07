%% Prep Sessions/Animals, prep variables for later blocks
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%

subDatT = import_xldat("D:\Data\Kelton\analyses\group_analyses","dat_include.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\Subiculum_RZ_Shift\rz_rand';
cd(groupSDir) 

cleanInds = subDatT.include ~= 1 | subDatT.sess_type ~= 3;
subDatT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;
sbase = 'subRwdShift_';
fname = [sbase 'data1'];

dbnsz = 0.05;
histoBnsz = 5;
binedges = 0:5:185;
binpos = 0.025:dbnsz:1.825;
wlen = 150;
r1pos = 0.1;    % 10 cm
r2pos = 1;      % 100cm
% vColors2 = [0.5 0.5 1; 0.75 0.75 1];
vColors2 = [.35 .35 .35; 1 .25 .25];

clear ps stats

%% Combine data
combine_rzShiftDat(subDatT,fname,groupSDir,'ca1',2); % sessType = 2

%% Load previously saved data
cd(groupSDir)
load(fname)

nMice = length(unique(recID(:,1)));
mID = unique(recID(:,1));
nTotal = size(recID,1);
r1posInd = find(binpos > r1pos,1);
r2posInd = find(binpos > r2pos,1);

vlFrstID = useCC & vlDat(:,1) <= 0.05;
vlLastID = useCC & vlDat(:,4) <= 0.05;
vlBothID = vlFrstID & vlLastID;
thFrstID = useCC & thDat(:,1) <= 0.05;
thLastID = useCC & thDat(:,4) <= 0.05;
thBothID = thFrstID & thLastID;
siFrstID = useCC & lcDat(:,1) <= 0.05;
siLastID = useCC & lcDat(:,5) <= 0.05;
siBothID = siFrstID & siLastID;
swrFrstID = useCC & rpDat(:,2) > 1;
swrLastID = useCC & rpDat(:,4) > 1;
swrBothID = swrFrstID & swrLastID;
bstFrstID = useCC & bsDat(:,1) > 0;
bstLastID = useCC & bsDat(:,2) > 0;
bstBothID = bstFrstID & bstLastID;

%% Behavior comparisons

uLckDI = [vertcat(bvDat.uPreLckDI), vertcat(bvDat.uPstLckDI)];
uRZVel = [vertcat(bvDat.uPreRZVel), vertcat(bvDat.uPstRZVel)];

[~,ps.bhv_PPLckDI,~,stats.bhv_PPLckDI] = ttest(uLckDI(:,1),uLckDI(:,2));
[~,ps.bhv_PPRZ1,~,stats.bhv_PPRZ1] = ttest(uRZVel(:,1),uRZVel(:,3));
[~,ps.bhv_PPRZ2,~,stats.bhv_PPRZ2] = ttest(uRZVel(:,2),uRZVel(:,4));

uLckDIFig = plot_barXmouse(uLckDI); ylim([-1 1.15])
text2bar(uLckDIFig,'Lick DI',ps.bhv_PPLckDI);

uRZ1VlFig = plot_barXmouse(uRZVel(:,[1,3])); ylim([0 50])
text2bar(uRZ1VlFig,'Velocity (cm/s) 30cm Pre RZ1',ps.bhv_PPRZ1);

uRZ2VlFig = plot_barXmouse(uRZVel(:,[2,4]));
text2bar(uRZ2VlFig,'Velocity (cm/s) 30cm Pre RZ2',ps.bhv_PPRZ2);

if saveFlag
    fsave(uLckDIFig,[sbase 'bhv_LckDI_bar'])
    fsave(uRZ1VlFig,[sbase 'bhv_RZ1Vl_bar'])
    fsave(uRZ2VlFig,[sbase 'bhv_RZ2Vl_bar'])
end

%% Behavior comparisons for last 20 laps in each RZ 
% uLckDI20 = [vertcat(bvDat.uPreLckDI20), vertcat(bvDat.uPstLckDI20)];
% uRZVel20 = [vertcat(bvDat.uPreRZVel20), vertcat(bvDat.uPstRZVel20)];
% 
% [~,ps.bhv_PPLckDI20,~,stats.bhv_PPLckDI20] = ttest(uLckDI20(:,1),uLckDI20(:,2));
% [~,ps.bhv_PPRZ120,~,stats.bhv_PPRZ120] = ttest(uRZVel20(:,1),uRZVel20(:,3));
% [~,ps.bhv_PPRZ220,~,stats.bhv_PPRZ220] = ttest(uRZVel20(:,2),uRZVel20(:,4));
% 
% uLckDIFig = plot_barXmouse(uLckDI20);
% ylabel('Lick DI ((RZ - AZ) / (RZ + AZ))'); ylim([-1 1])
% 
% uRZ1VlFig = plot_barXmouse(uRZVel20(:,[1,3]));
% ylabel('Velocity (cm/s) 30cm Pre RZ1')
% 
% uRZ2VlFig = plot_barXmouse(uRZVel20(:,[2,4]));
% ylabel('Velocity (cm/s) 30cm Pre RZ2')
% 
% if saveFlag
%     fsave(uLckDIFig,[sbase 'bhv_LckDI20_bar'])
%     fsave(uRZ1VlFig,[sbase 'bhv_RZ1Vl20_bar'])
%     fsave(uRZ2VlFig,[sbase 'bhv_RZ2Vl20_bar'])
% end

%% Averaged lick and velocity profiles

for i = 1:length(mID)
    uPreVel(i,:) = mean(bvDat(i).preVMap,'omitnan');
    uPstVel(i,:) = mean(bvDat(i).pstVMap,'omitnan');

    uPreLck(i,:) = mean(bvDat(i).preLMap,'omitnan');
    uPstLck(i,:) = mean(bvDat(i).pstLMap,'omitnan');
end

xcoords = [10 40; 100 130];
uVelF = plot_bhvrTraceCI(circshift(uPreVel,30,2),circshift(uPstVel,30,2),[0.45 0.45 0.45; 0 0 0],xcoords,[0 50]);
xticks([0 40 80 120 160]); xticklabels([-30 10 50 90 130]);
ylabel('Velocity (cm/s)'); ylim([0 50]); xlim([0 186])

uLckF = plot_bhvrTraceCI(circshift(uPreLck,10,2),circshift(uPstLck,10,2),[0.45 0.45 0.45; 0 0 0],xcoords,[0 12]);
xticks([0 40 80 120 160]); xticklabels([-30 10 50 90 130]);
ylabel('Lick rate (Hz)'); ylim([0 12]); xlim([0 186])

if saveFlag
    fsave(uVelF,[sbase 'bhv_meanVelPrePst'])
    fsave(uLckF,[sbase 'bhv_meanLckPrePst'])
end

%% Firing Rate stand vs run
% frDat: 1&3 = standing; 2&4 = running
nUseCC = sum(useCC);

[~,ps.fr_PPStnd_all,~,stats.fr_PPStnd_all] = ttest(frDat(useCC,1),frDat(useCC,3));
[~,ps.fr_PPRunn_all,~,stats.fr_PPRunn_all] = ttest(frDat(useCC,2),frDat(useCC,4));

bardat = [mean(frDat(useCC,1)), mean(frDat(useCC,2)); mean(frDat(useCC,3)), mean(frDat(useCC,4))];
xcoords = ones(sum(useCC),1);

frStndRunFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.15 0.27])
b1 = bar(bardat','FaceColor','flat');
b1(1).CData = vColors2(1,:);
b1(2).CData = vColors2(2,:);
errorbar([0.85 1.15 1.85 2.15],reshape(bardat,4,1),std(frDat)./sqrt(nUseCC),'k.')
% plot(xcoords-0.15,frDat(useCC,1),'ko')
% plot(xcoords+0.15,frDat(useCC,3),'ko')
% plot(xcoords+1-0.15,frDat(useCC,2),'ko')
% plot(xcoords+1+0.15,frDat(useCC,4),'ko')
xlim([0.5 2.5])
xticks(1:2); xticklabels({'Standing', 'Running'})
ylabel('Firing Rate (Hz)')
legend({'Familiar','Novel'},'Location','northwest')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    fsave(frStndRunFig,[sbase 'fr_bar'])
end

%% Burst Rate pre/post
% bsDat: 1 = familiar; 2 = novel

[~,ps.bs_BI_both,~,stats.bs_BI_both] = ttest(bsDat(bstBothID,1), bsDat(bstBothID,2));

% Units with bursts in both epochs
bsFig = plotBar2(bsDat(bstBothID,1),bsDat(bstBothID,2)); ylim([0 6])
text2bar(bsFig,'Burst Index',ps.bs_BI_both);
bsFigZoom = plotBar2(bsDat(bstBothID,1),bsDat(bstBothID,2)); ylim([0 1.2])

if saveFlag
    fsave(bsFig,[sbase 'bs_BI_both'])
    fsave(bsFigZoom,[sbase 'bs_BI_both_zoom'])
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

%% Spatial Information, peak rate
% lcDat: 1&5 = sig.; 2&6 = SI; 3&7 = pkRate; 4&8 = pkLoc

[~,ps.lc_PPSI_both,~,stats.lc_PPSI_both] = ttest(lcDat(siBothID,2), lcDat(siBothID,6));
[~,ps.lc_PPPk_both,~,stats.lc_PPPk_both] = ttest(lcDat(siBothID,3), lcDat(siBothID,7));
[~,ps.lc_PPSI_eith,~,stats.lc_PPSI_eith] = ttest2(lcDat(siFrstID & ~siLastID,2),lcDat(siLastID & ~siFrstID,6));
[~,ps.lc_PPPk_eith,~,stats.lc_PPPk_eith] = ttest2(lcDat(siFrstID & ~siLastID,3),lcDat(siLastID & ~siFrstID,7));

% Summary of quantities
[siPie, siprcts] = prepostPie(siFrstID,siLastID,useCC);
title("Sig. Spatial Info. units");
siBothCounts = [groupcounts(recID(siFrstID,1),[mID;mID(end)+1],'IncludeEmptyGroups',true) groupcounts(recID(siLastID,1),[mID;mID(end)+1],'IncludeEmptyGroups',true)];
siBothRatio = siBothCounts ./ [groupcounts(recID(useCC,1)) groupcounts(recID(useCC,1))];
[~,ps.lc_SICt_both,~,stats.lc_SICt_both] = ttest(siBothCounts(:,1),siBothCounts(:,2));
lcBothCtFig = plot_barXmouse(siBothRatio);
text2bar(lcBothCtFig,'P(Sig. Spatial Info.)',ps.lc_SICt_both)

% Units modulated in both phases
lcBothSIFig = plotBar2(lcDat(siBothID,2),lcDat(siBothID,6)); ylim([0 4])
text2bar(lcBothSIFig,'Spatial Information (Bits/spike)',ps.lc_PPSI_both);
lcBothPkFig = plotBar2(lcDat(siBothID,3),lcDat(siBothID,7));
text2bar(lcBothPkFig,'Peak Field FR (Hz)',ps.lc_PPPk_both);

% For units only modulated in one task phase
lcEithSIFig = plotBar2(lcDat(siFrstID & ~siLastID,2),lcDat(siLastID & ~siFrstID,6));
text2bar(lcEithSIFig,'Spatial Information (Bits/spike)',ps.lc_PPSI_eith); xticklabels({'F-Only', 'N-Only'})
lcEithPkFig = plotBar2(lcDat(siFrstID & ~siLastID,3),lcDat(siLastID & ~siFrstID,7));
text2bar(lcEithPkFig,'Peak Field FR (Hz)',ps.lc_PPPk_eith); xticklabels({'F-Only', 'N-Only'})

if saveFlag
    fsave(siPie,[sbase 'si_Mod_pie'])
    fsave(lcBothSIFig,[sbase 'si_both_bar'])
    fsave(lcBothPkFig,[sbase 'pk_both_bar'])
    fsave(lcEithSIFig,[sbase 'si_eith_bar'])
    fsave(lcEithPkFig,[sbase 'pk_eith_bar'])
    fsave(lcBothCtFig,[sbase 'si_Mod_both_bar'])
end

%% Waterfall by position
binedges = 0:0.05:1.85;
nBins = length(binedges) - 1;

[spBothPreSortPreFig,tmpMap,sortPre] = plot_unitWaterfall(circshift(lcMap(siBothID,1:length(binpos)),18,2),binedges,0,1,0);
xticks([1, round(nBins)/2, nBins]); xticklabels([-90, 0, 90])
plot([r2posInd r2posInd],[0 sum(siBothID)],'k--','LineWidth',2); title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
[spBothPreSortPreHisto,pkMapPre] = plot_unitPkHisto(tmpMap,binedges*100,1);
xticks([binedges(1),binedges(round(nBins/2)),binedges(end)]*100); xticklabels([-90, 0, 90])
plot([r2pos r2pos]*100,[0 0.11],'k--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.1])
[ps.lc_bothPrePkUniformity, stats.lc_bothPrePkUniformity] = pkChi2(pkMapPre,binedges);
text2bar(spBothPreSortPreHisto,'',ps.lc_bothPrePkUniformity,0.9)

[spBothPstSortPstFig,tmpMap] = plot_unitWaterfall(lcMap(siBothID,length(binpos)+1:end),binedges,0,1,0);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2); title('Novel RZ, sort Novel'); xlabel('Track Position (cm)')
[spBothPstSortPstHisto,pkMapPst] = plot_unitPkHisto(tmpMap,binedges*100,1);
plot([r2pos r2pos]*100,[0 0.11],'r--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.1])
[ps.lc_bothPstPkUniformity, stats.lc_bothPstPkUniformity] = pkChi2(pkMapPst,binedges);
text2bar(spBothPstSortPstHisto,'',ps.lc_bothPstPkUniformity,0.9)

[spBothPstSortPreFig,tmpMap] = plot_unitWaterfall(lcMap(siBothID,length(binpos)+1:end),binedges,sortPre,1,0);
plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2); title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')
spBothPstSortPreHisto = plot_unitPkHisto(tmpMap,binedges*100,1);
plot([r2pos r2pos]*100,[0 0.11],'r--','LineWidth',2); xlabel('Track Position (cm)'); ylim([0 0.1])

siBothPkDistroHisto = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.14])
plot([binedges(1:end-1) + 0.5*diff(binedges(1:2))]*100,sum(pkMapPre)./sum(pkMapPre,'all'),'Color',vColors2(1,:),'LineWidth',2)
plot([binedges(1:end-1) + 0.5*diff(binedges(1:2))]*100,sum(pkMapPst)./sum(pkMapPst,'all'),'Color',vColors2(2,:),'LineWidth',2)
set(gca,'Position',[0.11 0.17 0.8 0.80]); xlim([0 186]); ylim([0 0.1])
xticks(100*[binedges(1), binedges(round(nBins/2)), binedges(nBins)]); 
set(gca,'FontSize',12,'FontName','Arial')

% [spEithPreSortPreFig,~,sortPre] = plot_unitWaterfall(lcMap(siFrstID & ~siLastID,1:length(binpos)),binedges,0,1,0);
% plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Familiar RZ, sort Familiar'); xlabel('Track Position (cm)')
% 
% [spEithPstSortPreFig] = plot_unitWaterfall(lcMap(siFrstID & ~siLastID,length(binpos)+1:end),binedges,sortPre);
% plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Novel RZ, sort Familiar'); xlabel('Track Position (cm)')

% [spEithPstSortPstFig,tmpMap,Sortpst] = plot_unitWaterfall(lcMap(siLastID & ~siFrstID,length(binpos)+1:end),binedges);
% plot([r2posInd r2posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Novel RZ, sort Novel'); xlabel('Track Position (cm)')
% [spEithPstSortPstHisto,tmpPks] = plot_unitPkHisto(tmpMap,binedges*100);

% probs = sum(tmpPks) ./ sum(tmpPks,'all');
% expec = size(tmpPks,1) ./ size(tmpPks,2) * ones(size(tmpPks,2),1);
% [~,p,stats] = chi2gof(cumsum(probs),'Expected',expec,'Edges',linspace(0,1,length(binedges)))
% xlabel('Track Position (cm)'); ylim([0 0.11])

% [spEithPreSortPstFig] = plot_unitWaterfall(lcMap(siLastID & ~siFrstID,1:length(binpos)),binedges,Sortpst);
% plot([r1posInd r1posInd],[0 sum(siBothID)],'r--','LineWidth',2)
% title('Familiar RZ, sort Novel'); xlabel('Track Position (cm)')

% wfCBarF = plotColorbar([0 1],'parula');

if saveFlag
    fsave(spBothPreSortPreFig,[sbase 'sp_both_pre_sortPre'])
    fsave(spBothPstSortPstFig,[sbase 'sp_both_pst_sortPst'])
    fsave(spBothPstSortPreFig,[sbase 'sp_both_pst_sortPre'])
    fsave(spBothPreSortPreHisto,[sbase 'sp_both_pre_sortPre_hist'])
    fsave(spBothPstSortPstHisto,[sbase 'sp_both_pst_sortPst_hist'])
    fsave(spBothPstSortPreHisto,[sbase 'sp_both_pst_sortPre_hist'])
    fsave(siBothPkDistroHisto,[sbase 'sp_both_hist'])
    % saveas(spEithPreSortPreFig,[sbase 'sp_eith_pre_sortPre'],'png')
    % saveas(spEithPstSortPreFig,[sbase 'sp_eith_pst_sortPre'],'png')
    % saveas(spEithPstSortPstFig,[sbase 'sp_eith_pst_sortPst'],'png')
    % saveas(spEithPreSortPstFig,[sbase 'sp_eith_pre_sortPst'],'png')
end

%% Bar summary of SI, Vel, and Theta tuning

bardat = [siprcts; thprcts; vlprcts; swprcts];

figure; hold on; 
set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.39])
b = barh(bardat,'stacked');
b(1).FaceColor = [0.75 0.75 1];
b(2).FaceColor = vColors2(1,:);
b(3).FaceColor = vColors2(2,:);
b(4).FaceColor = [1 1 1];
yticks(1:4); yticklabels({'Spatial','Theta','Velocity','SPW-R'});
xlabel('Percentage Modulated'); ylim([0.5 4.5])
legend({'Both','Fam. only','Nov. only','Neither'})
box off
set(gca,'FontSize',16,'FontName','Arial')

%% SI over 10 trial blocks
% 
% cmapcool = cool(nMice);
% siEpochFig = figure; hold on
% for i = 1:nMice
%     subID = siBothID(recID(:,1) == mID(i));
%     siStr(i).meanPreSI = mean(siStr(i).preBlockSI(subID,:));
%     nEpochsPre(i) = size(siStr(i).preBlockSI,2);
%     plot(1:nEpochsPre(i),siStr(i).meanPreSI,'-o','Color',cmapcool(i,:))
%     nEpochsPst(i) = size(siStr(i).pstBlockSI,2);
% end
% maxEpoch = max(nEpochsPre);
% uPreMat = NaN(size(recID,1),maxEpoch);
% uPstMat = NaN(size(recID,1),max(nEpochsPst));
% 
% ct = 1;
% for i = 1:nMice
%     nUnits = sum(recID(:,1) == mID(i));
%     uPreMat(ct:nUnits+ct-1,1:nEpochsPre(i)) = siStr(i).preBlockSI;
%     uPstMat(ct:nUnits+ct-1,1:nEpochsPst(i)) = siStr(i).pstBlockSI;
%     ct = ct+nUnits;
% end
% for i = 1:nMice
%     subID = siBothID(recID(:,1) == mID(i));
%     siStr(i).meanPstSI = mean(siStr(i).pstBlockSI(subID,:));
%     nPst(i) = size(siStr(i).pstBlockSI,2);
%     plot(maxEpoch+1:nPst(i)+maxEpoch,siStr(i).meanPstSI,'-o','Color',cmapcool(i,:))
% end
% 
% errorbar(1:maxEpoch,mean(uPreMat(siBothID,:),'omitnan'),std(uPreMat(siBothID,:),'omitnan')./sqrt(nTotal),'k-','LineWidth',2)
% errorbar(maxEpoch+1:max(nEpochsPst)+maxEpoch,mean(uPstMat(siBothID,:),'omitnan'),std(uPstMat(siBothID,:),'omitnan')./sqrt(nTotal),'k-','LineWidth',2)
% plot([maxEpoch+0.5 maxEpoch+0.5],[0 max(ylim)],'k--')
% xlabel('10-Trial block #')
% ylabel('Spatial Info. (bits/spike)')
% % ylim([0 1])
% set(gca,'FontSize',12,'FontName','Arial')
% 
% if saveFlag
%     saveas(siEpochFig,[sbase 'si_both_subEpoch'],'png')
% end

%% Field peak spatial distribution

% % Familiar vs Novel RZ correlation
% xrand = 2*rand(size(lcDat(siBothID,4)))-1;
% yrand = 2*rand(size(lcDat(siBothID,4)))-1;
% pkPosFig = figure; axis square; hold on
% plot(100*binpos(lcDat(siBothID,4))+xrand',100*binpos(lcDat(siBothID,8))+yrand','k.','MarkerSize',10)
% plot([0 100*binpos(end)],[0 100*binpos(end)],'k--')
% plot([r1pos r1pos]*100,[0 100*binpos(end)],'r--',[0 100*binpos(end)],[r2pos r2pos]*100,'r--')
% xlabel('Absolute Peak Loc. (cm) Familiar'); xlim([0 100*binpos(end)])
% ylabel('Absolute Peak Loc. (cm) Novel');    ylim([0 100*binpos(end)])
% set(gca,'FontSize',12,'FontName','Arial')
% 
% % Familiar vs Novel RZ correlation heatmap
% fieldDstPrePst = histcounts2(binpos(lcDat(siBothID,8)),binpos(lcDat(siBothID,4)),binedges,binedges);
% pkPosMapFig = figure; hold on; axis square
% imagesc(fieldDstPrePst)
% colormap("parula")
% cbar = colorbar; %clim([prctile(fieldDstPrePst,1,'all'), prctile(fieldDstPrePst,99,'all')]);
% plot([0 length(binpos)],[0 length(binpos)],'w--')
% plot([r1posInd r1posInd],[0 length(binpos)],'r--',[0 length(binpos)],[r2posInd r2posInd],'r--')
% xlabel('Absolute Peak Loc. (cm) Familiar'); ylabel('Absolute Peak Loc. (cm) Novel')
% ylabel(cbar,'Count','FontSize',12,'Rotation',90); 
% xlim([0 length(binedges)]); ylim([0 length(binedges)])
% xticks(1:10:length(binpos)); yticks(1:10:length(binpos))
% xticklabels(binedges(1:10:end)*100); yticklabels(binedges(1:10:end)*100)
% set(gca,'FontSize',12,'FontName','Arial','YDir','normal','XDir','normal')
% 
% if saveFlag
%     saveas(pkPosFig,[sbase 'lc_PeakComp'],'png')
%     saveas(fieldDstPrePst,[sbase 'lc_PeakMap'],'png')
% end

%% Finding Track Relative or Reward Relative cells

spPk1 = binpos(lcDat(:,4));
spPk2 = binpos(lcDat(:,8));

dth = 0.3;  % Distance from peak threshold (m)
drz = r2pos - r1pos;
trCells = (spPk2 - spPk1 < dth & spPk2 - spPk1 > -dth)' & siBothID; % threshold 30cm
rrCells = ((spPk2 - spPk1 > drz-dth & spPk2 - spPk1 < drz+dth)' & siBothID) | ((spPk2 - spPk1 < -drz+dth & spPk2 - spPk1 > -drz-dth)' & siBothID);
irCells = siBothID & ~trCells & ~rrCells;
xrand = 2*rand(size(lcDat(siBothID,4)))-1;
yrand = 2*rand(size(lcDat(siBothID,4)))-1;

pkPosPatchFig = figure; hold on; axis square
set(gcf,'units','normalized','position',[0.4 0.35 0.24 0.39])
patch(100*[0 dth 1.85 1.85 1.85-dth 0 0],100*[0 0 1.85-dth 1.85 1.85 dth 0],...
    'b','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
patch(100*[0 1.25 0.65 0 0],100*[0.6 1.85 1.85 1.2 0.6],...
    'r','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
patch(100*[1.2 1.85 1.85 0.6 1.2],100*[0 0.65 1.25 0 0],...
    'r','FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
plot(100*spPk1(siBothID)+xrand',100*spPk2(siBothID)+yrand','k.','MarkerSize',10)
plot([0 100*binpos(end)],[0 100*binpos(end)],'k--')
plot([r1pos r1pos]*100,[0 100*binpos(end)],'r--',[0 100*binpos(end)],[r2pos r2pos]*100,'r--')
xlabel('Absolute Peak Loc. (cm) Familiar'); xlim([0 100*binpos(end)])
ylabel('Absolute Peak Loc. (cm) Novel');    ylim([0 100*binpos(end)])
set(gca,'FontSize',16,'FontName','Arial')

nUseCC = sum(useCC);

nBoth = sum(siBothID);
nPre = sum(siFrstID);
nPost = sum(siLastID);
nNot = sum(not(siFrstID | siLastID) & useCC);

prcts = [sum(trCells), sum(irCells), sum(rrCells), nPre - nBoth, nPost - nBoth, nNot] ./ sum(useCC); 

cMap = [.75, .75, 1; .75, .5, .75; 1, .75, .75; .35, .35, .35; 1, .25, .25; 1, 1, 1];
trirrrPie = figure;
p = piechart(prcts*100,["TR-cells","IR-cells","RR-cells","Familiar-only","Novel-only","Neither"]);
p.LabelStyle = 'namedata';
colororder(cMap)

if saveFlag
    fsave(pkPosPatchFig,[sbase 'lc_PeakComp_patch'])
    fsave(trirrrPie,[sbase 'lc_PeakComp_pie'])
end

%% Distributions over linearized space
binedges = 0:0.05:1.85;

activeCells = rrCells; % siBothID; rrCells; trCells
fieldDstPre = histcounts(spPk1(activeCells),binedges);
fieldDstPst = histcounts(spPk2(activeCells),binedges);

prepstFieldDstFig = plotDistroHisto(fieldDstPre,fieldDstPst,binpos,[r1pos r2pos]);
xlabel('Track Position (cm)')

% Aligns peaks with 90cm as reward
trackLen = max(binedges);
shiftR1 = trackLen/2 - r1pos;
fieldAlignR1 = mod(binpos(lcDat(activeCells,4))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
shiftR2 = trackLen/2 - r2pos;
fieldAlignR2 = mod(binpos(lcDat(activeCells,8))+shiftR2,trackLen)+0.5*dbnsz;
shiftbins = binedges - trackLen/2;

fieldDstPre_Align = histcounts(fieldAlignR1,binedges);
fieldDstPst_Align = histcounts(fieldAlignR2,binedges);

rzDstFig = plotDistroHisto(fieldDstPre_Align,fieldDstPst_Align,binpos-trackLen/2,0);
xlabel('Distance to RZ (cm)')

deltaField_RZAlign = fieldAlignR2 - fieldAlignR1;
uFieldAlign = mean([fieldAlignR1; fieldAlignR2],1);
uFieldAlignDistro = histcounts(uFieldAlign-trackLen/2,shiftbins);
circAlignNeg = deltaField_RZAlign < -trackLen/2;
circAlignPos = deltaField_RZAlign > trackLen/2;
deltaField_RZAlign(circAlignNeg) = -(deltaField_RZAlign(circAlignNeg) + trackLen);     % When new field back-shifts
deltaField_RZAlign(circAlignPos) = -(deltaField_RZAlign(circAlignPos) - trackLen);     % When new field forward-shifts

[~,ps.lc_DltLc_both,~,stats.lc_DltLc_both] = ttest(deltaField_RZAlign);

deltaField_Distro = histcounts(deltaField_RZAlign,shiftbins);

fieldShiftRZFig = figure; hold on
bar((shiftbins(1:end-1)+0.5*dbnsz)*100,deltaField_Distro/sum(deltaField_Distro),'FaceColor',[0.25 0.15 1],'HandleVisibility','off')
plot([0 0]*100,[0 max(deltaField_Distro/sum(deltaField_Distro),[],'all')+0.02],'--','Color',[.5 .5 .5],'HandleVisibility','off')
plot(mean(deltaField_RZAlign),max(deltaField_Distro/sum(deltaField_Distro),[],'all')+0.02,'v','Color',[0.25 0.15 1])
xlabel('\Delta RZ-aligned Novel - Familiar (cm)'); ylabel('Probability')
xlim([-trackLen/2-dbnsz trackLen/2+dbnsz]*100)
legend({['mean = ' num2str(mean(deltaField_RZAlign))]})
set(gca,'FontSize',12,'FontName','Arial')

%% Shuffle novel RZ peak locations and calculate global confidence bands 
clear deltaField_DistroJit
for i = 1:250
    rShift = 1 + randi(length(binedges) - 2,sum(activeCells),1);
    jitPost = mod(lcDat(activeCells,8) + rShift, length(binedges)-1)+1;
    fieldAlignJit = mod(binpos(jitPost)+shiftR2,trackLen)+0.5*dbnsz;
    fieldDstJit_Align = histcounts(fieldAlignJit,binedges);
    deltaField_RZJit = fieldAlignJit - fieldAlignR1;
    circAlignNeg = deltaField_RZJit < -trackLen/2;
    circAlignPos = deltaField_RZJit > trackLen/2;
    deltaField_RZJit(circAlignNeg) = -(deltaField_RZJit(circAlignNeg) + trackLen);     % When new field back-shifts
    deltaField_RZJit(circAlignPos) = -(deltaField_RZJit(circAlignPos) - trackLen);     % When new field forward-shifts
    deltaField_DistroJit(:,i) = histcounts(deltaField_RZJit,shiftbins);

    % uFieldAlignJit = mean([fieldAlignR1; fieldAlignJit],1);
    % uFieldDistroJit(:,i) = histcounts(uFieldAlignJit - trackLen/2,shiftbins);
    uFieldDistroJit(:,i) = histcounts(fieldAlignJit-trackLen/2,shiftbins);
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

activeUnits = siBothID;     % siBothID (allSICells); trCells; rrCells;
nBins = length(binpos);

posNormPre = lcMap(activeUnits,1:nBins) ./ max(lcMap(activeUnits,1:nBins),[],2);
posNormPst = lcMap(activeUnits,nBins+1:end) ./ max(lcMap(activeUnits,nBins+1:end),[],2);

pvPrePst = corr(posNormPst,posNormPre);
pvPrePstF = plot_pvcorr(pvPrePst);

% Compare identity to off-diagonal
idMat = logical(eye(size(pvPrePst)));
dgMat = logical(spdiags([1 1],[-round(nBins/2) round(nBins/2)],nBins,nBins));

allPreEvn = [];
allPreOdd = [];
allPstEvn = [];
allPstOdd = [];

for i = 1:nMice
    tmpUnits = activeUnits & recID(:,1) == mID(i);

    posNormPre = lcMap(tmpUnits,1:nBins) ./ max(lcMap(tmpUnits,1:nBins),[],2);
    posNormPst = lcMap(tmpUnits,nBins+1:end) ./ max(lcMap(tmpUnits,nBins+1:end),[],2);
    try
        pvPrePst = corr(posNormPst,posNormPre);
    catch
        pvPrePst = NaN(nBins,nBins);
    end
    uPVCorrID(i) = mean(pvPrePst(idMat),'all');
    uPVCorrDG(i) = mean(pvPrePst(dgMat),'all');

    allPreEvn = [allPreEvn; pvStr(i).preEvn];
    allPreOdd = [allPreOdd; pvStr(i).preOdd];
    allPstEvn = [allPstEvn; pvStr(i).pstEvn];
    allPstOdd = [allPstOdd; pvStr(i).pstOdd];

    pvPreOddEvn = corr(pvStr(i).preOdd, pvStr(i).preEvn);
    pvPstOddEvn = corr(pvStr(i).pstOdd, pvStr(i).pstEvn);

    uPVCPreOddEvnID(i) = mean(pvPreOddEvn(idMat),'all','omitnan');
    uPVCPreOddEvnDG(i) = mean(pvPreOddEvn(dgMat),'all','omitnan');
    uPVCPstOddEvnID(i) = mean(pvPstOddEvn(idMat),'all','omitnan');
    uPVCPstOddEvnDG(i) = mean(pvPstOddEvn(dgMat),'all','omitnan');
end

pvPreOddEvn = corr(allPreOdd,allPreEvn,'rows','complete');
pvPreOddEvnF = plot_pvcorr(pvPreOddEvn);
xlabel("Position (cm) even laps"); ylabel("Position (cm) odd laps")

pvPstOddEvn = corr(allPstOdd,allPstEvn,'rows','complete');
pvPstOddEvnF = plot_pvcorr(pvPstOddEvn);
xlabel("Position (cm) even laps"); ylabel("Position (cm) odd laps")

[~,ps.lc_pv_PrePst_idVdg,~,stats.lc_pv_PrePst_idVdg] = ttest(uPVCorrID,uPVCorrDG);
[~,ps.lc_pv_PreOddEvn_idVdg,~,stats.lc_pv_PreOddEvn_idVdg] = ttest(uPVCPreOddEvnID,uPVCPreOddEvnDG);
[~,ps.lc_pv_PstOddEvn_idVdg,~,stats.lc_pv_PstOddEvn_idVdg] = ttest(uPVCPstOddEvnID,uPVCPstOddEvnDG);

pvPrePstCompF    = plot_PVCorrComp(uPVCorrID, uPVCorrDG, ps.lc_pv_PrePst_idVdg);
pvPreOddEvnCompF = plot_PVCorrComp(uPVCPreOddEvnID, uPVCPreOddEvnDG, ps.lc_pv_PreOddEvn_idVdg);
pvPstOddEvnCompF = plot_PVCorrComp(uPVCPstOddEvnID, uPVCPstOddEvnDG, ps.lc_pv_PstOddEvn_idVdg);

if saveFlag
    fsave(pvPrePstF,[sbase 'lc_pv_corr_allSICells_prepost'])
    fsave(pvPreOddEvnF,[sbase 'lc_pv_corr_allSICells_preoddeven'])
    fsave(pvPstOddEvnF,[sbase 'lc_pv_corr_allSICells_pstoddeven'])
    fsave(pvPrePstCompF,[sbase 'lc_pv_comp_allSICells_prepost'])
    fsave(pvPreOddEvnCompF,[sbase 'lc_pv_comp_allSICells_preoddeven'])
    fsave(pvPstOddEvnCompF,[sbase 'lc_pv_comp_allSICells_pstoddeven'])
end

%% PV across time

lapcutoff = 65;
for i = 1:nMice
    nTrialsPre(i) = size(pvStr(i).preBlockPV,1);
    nTrialsPst(i) = size(pvStr(i).pstBlockPV,1);
end
maxTr = [max(nTrialsPre) max(nTrialsPst)];
pvXlapPrePre = NaN(maxTr(1),nMice);
pvXlapPrePst = NaN(maxTr(1),nMice);
pvXlapPstPre = NaN(maxTr(2),nMice);
pvXlapPstPst = NaN(maxTr(2),nMice);

for i = 1:nMice
    pvXlapPrePre(1:nTrialsPre(i),i) = pvStr(i).preBlockPV(:,1);
    pvXlapPrePst(1:nTrialsPre(i),i) = pvStr(i).preBlockPV(:,2);
    pvXlapPstPre(1:nTrialsPst(i),i) = pvStr(i).pstBlockPV(:,1);
    pvXlapPstPst(1:nTrialsPst(i),i) = pvStr(i).pstBlockPV(:,2);
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
plot_CIs(1:lapcutoff,preCIup,preCIdn,vColors2(1,:)); % 1:maxTr(1)
plot_CIs(1:lapcutoff,pstCIup,pstCIdn,vColors2(2,:)); % 1:maxTr(2)
plot(1:lapcutoff,nanmean(pvXlapPrePre'),'Color',vColors2(1,:),'LineWidth',2)
plot(1:lapcutoff,nanmean(pvXlapPstPst'),'Color',vColors2(2,:),'LineWidth',2)
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
plot_CIs(1:lapcutoff,preCIup,preCIdn,vColors2(1,:))
plot_CIs(1:lapcutoff,pstCIup,pstCIdn,vColors2(2,:))
plot(1:lapcutoff,nanmean(vcXlapPre'),'Color',vColors2(1,:),'LineWidth',2)
plot(1:lapcutoff,nanmean(vcXlapPst'),'Color',vColors2(2,:),'LineWidth',2)
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
plot_CIs(1:lapcutoff,preCIup,preCIdn,vColors2(1,:))
plot_CIs(1:lapcutoff,pstCIup,pstCIdn,vColors2(2,:))
plot(1:lapcutoff,nanmean(lDIXlapPre'),'Color',vColors2(1,:),'LineWidth',2)
plot(1:lapcutoff,nanmean(lDIXlapPst'),'Color',vColors2(2,:),'LineWidth',2)
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
plot(vcXpvcPrePre(:,1),vcXpvcPrePre(:,2),'o','Color',vColors2(1,:))
plot(vcXpvcPstPst(:,1),vcXpvcPstPst(:,2),'o','Color',vColors2(2,:))
plot(vcXpvcPrePre(:,1),mdlVPreXPVCprepre.ypred,'Color',vColors2(1,:),'LineWidth',2)
plot(vcXpvcPstPst(:,1),mdlVPreXPVCpstpst.ypred,'Color',vColors2(2,:),'LineWidth',2)
xlim([-1 1]); xlabel('Velocity Correlation');
ylim([0 1]); ylabel('PV Correlation');
legend({'F-F Vel -> F-F PVC','N-N Vel -> N-N PVC'},'location','ne')
set(gca,'FontSize',16,'FontName','Arial')
xlims = xlim;
ylims = ylim;
text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims),  ['R = ' num2str(mdlVPreXPVCprepre.r, 3)], 'Color', vColors2(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlVPreXPVCprepre.p, 3)], 'Color', vColors2(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.25*diff(ylims), ['R = ' num2str(mdlVPreXPVCpstpst.r, 3)], 'Color', vColors2(2,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.3*diff(ylims),  ['p = ' num2str(mdlVPreXPVCpstpst.p, 3)], 'Color', vColors2(2,:), 'FontSize', 12)

if saveFlag
    fsave(vcXpvcF,[sbase '_vCorrXPVC'])
end

%% SPWR Data
% rpDat: 2&4 = sig.; 1&3 = participation; 
% rpRat: 1&2 = rate by mouse

% Summary of quantities
[tmpSWPie, swprcts] = prepostPie(swrFrstID,swrLastID,useCC);
title("Sig. SPW-R Mod. units");
frstGrpAdd = [mID; recID(swrFrstID,1)]; % Need to add all mice then subtract 1 from groupCounts otherwise it is blind to missing mice
lastGrpAdd = [mID; recID(swrLastID,1)];
swrBothCounts = [groupcounts(firstGrpAdd)-1 groupcounts(lastGrpAdd)-1];
swrBothRatio = swrBothCounts ./ [groupcounts(recID(useCC,1)) groupcounts(recID(useCC,1))];
[~,ps.swr_ModCt_both,~,stats.swr_ModCt_both] = ttest(swrBothCounts(:,1),swrBothCounts(:,2));
swrModCtFig = plot_barXmouse(swrBothRatio);
text2bar(swrModCtFig,'P(Sig. SWR Mod)',ps.swr_ModCt_both);

[~,ps.swr_PPR_both,~,stats.swr_PPR_both] = ttest(rpRat(:,1),rpRat(:,2));
[~,ps.swr_PPP_both,~,stats.swr_PPP_both] = ttest(rpDat(swrBothID,1),rpDat(swrBothID,3));
[~,ps.swr_PPP_eith,~,stats.swr_PPP_eith] = ttest2(rpDat(swrFrstID & ~swrLastID,1),rpDat(swrLastID & ~swrFrstID,3));
[~,ps.swr_PPPre_trrr,~,stats.swr_PPPre_trrr] = ttest2(rpDat(trCells,1),rpDat(rrCells,1));
[~,ps.swr_PPPst_trrr,~,stats.swr_PPPst_trrr] = ttest2(rpDat(trCells,3),rpDat(rrCells,3));

rpRatFig = plot_barXmouse(rpRat);
text2bar(rpRatFig,'SPW-R Rate (Hz)',ps.swr_PPR_both);

% Units modulated in both phases
swPartcpFig = plotBar2(rpDat(swrBothID,1),rpDat(swrBothID,3));
text2bar(swPartcpFig,'P(SPW-R Participation)',ps.swr_PPP_both);

% For units only modulated in one task phase
swPartcpEithFig = plotBar2(rpDat(swrFrstID & ~swrLastID,1),rpDat(swrLastID & ~swrFrstID,3));
text2bar(swPartcpEithFig,'P(SPW-R Participation)',ps.swr_PPP_both);

if saveFlag
    fsave(tmpSWPie,[sbase 'swr_Mod_pie'])
    fsave(rpRatFig,[sbase 'swr_Rate_bar'])
    fsave(swPartcpFig,[sbase 'swr_Partcp_bar'])
    fsave(swrModCtFig,[sbase 'swr_ModCt_bar'])
    fsave(swPartcpEithFig,[sbase 'swr_Partcp_eith'])
end

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
plot(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(pkMapPre)./sum(pkMapPre,'all'),'Color',vColors2(1,:))
plot(binedges(1:end-1) + 0.5*diff(binedges(1:2)),sum(pkMapPst)./sum(pkMapPst,'all'),'Color',vColors2(2,:))
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

mdlLPreXspwr = get_linfit(uLckDI(:,1),swrBothRatio(:,1));
mdlLPstXspwr = get_linfit(uLckDI(:,2),swrBothRatio(:,2));

lckXspwrF = figure; hold on
plot(uLckDI(:,1),swrBothRatio(:,1),'o','Color',vColors2(1,:))
plot(uLckDI(:,2),swrBothRatio(:,2),'o','Color',vColors2(2,:))
plot(uLckDI(:,1),mdlLPreXspwr.ypred,'Color',vColors2(1,:),'LineWidth',2)
plot(uLckDI(:,2),mdlLPstXspwr.ypred,'Color',vColors2(2,:),'LineWidth',2)
xlim([-1 1]); xlabel('Lick DI');
ylim([0 1]); ylabel('P(SPWR-mod)');
legend({'Familiar','Novel'},'location','ne')
set(gca,'FontSize',16,'FontName','Arial')
xlims = xlim;
ylims = ylim;
text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims),  ['R = ' num2str(mdlLPreXspwr.r, 3)], 'Color', vColors2(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlLPreXspwr.p, 3)], 'Color', vColors2(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.25*diff(ylims), ['R = ' num2str(mdlLPstXspwr.r, 3)], 'Color', vColors2(2,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.3*diff(ylims),  ['p = ' num2str(mdlLPstXspwr.p, 3)], 'Color', vColors2(2,:), 'FontSize', 12)

% Compare delta of lckDI and delta of P(SPWR-Mod)
mdlDeltaLPreXspwr = get_linfit(uLckDI(:,2)-uLckDI(:,1),swrBothRatio(:,2)-swrBothRatio(:,1));
dlckXdspwrF = figure; hold on
plot(uLckDI(:,2)-uLckDI(:,1),swrBothRatio(:,2)-swrBothRatio(:,1),'ko')
plot(uLckDI(:,2)-uLckDI(:,1),mdlDeltaLPreXspwr.ypred,'Color','k','LineWidth',2)
xlabel('\Delta Lick DI');
ylabel('\Delta P(SPWR-mod)');
set(gca,'FontSize',16,'FontName','Arial')
xlims = xlim;ylims = ylim;
text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims),  ['R = ' num2str(mdlDeltaLPreXspwr.r, 3)], 'Color', vColors2(1,:), 'FontSize', 12)
text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlDeltaLPreXspwr.p, 3)], 'Color', vColors2(1,:), 'FontSize', 12)

if saveFlag
    fsave(lckXspwrF,[sbase 'swrXlck'])
    fsave(dlckXdspwrF,[sbase 'swrDXlckD'])
end

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
for i = 1:length(mID)
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

function [fhandle] = plot_bhvrTraceCI(dat1,dat2,vcolors,xcoords,ycoords)

[ciup1, cidn1] = get_CI(dat1);
[ciup2, cidn2] = get_CI(dat2);

bnpos = linspace(0,185,size(dat1,2));

xcoordsPre = xcoords(1,:);
xcoordsPst = xcoords(2,:);

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.5 0.22 0.23])

patch([xcoordsPre,fliplr(xcoordsPre)],[ycoords(1) ycoords(1) ycoords(2) ycoords(2)],[0 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
patch([xcoordsPst,fliplr(xcoordsPst)],[ycoords(1) ycoords(1) ycoords(2) ycoords(2)],[1 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
plot([xcoordsPre(2) xcoordsPre(2)],ycoords,'--','Color',vcolors(2,:)); plot([xcoordsPst(2) xcoordsPst(2)],ycoords,'--','Color',vcolors(2,:))

plot_CIs(bnpos,ciup1,cidn1,[.8 .8 .8])
plot(bnpos,mean(dat1),'Color',vcolors(2,:),'LineWidth',2)
plot_CIs(bnpos,ciup2,cidn2,[1 .8 .8])
plot(bnpos,mean(dat2),'Color',[1 0 0],'LineWidth',2)
xlabel('Position (cm)')
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
set(gcf,'units','normalized','position',[0.4 0.35 0.1 0.2])
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

