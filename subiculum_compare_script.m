% Compare and summarize ephys findings across subiculum sessions

%% Designate Analysis save location
sdir = 'D:\Data\Kelton\analyses\group_analyses\Subiculum_KW006-8';

%% Load in data
load('D:\Data\Kelton\analyses\KW006\KW006_07302024_rec_D1_Sub\KW006_07302024_rec_D1_Sub_root')
load('D:\Data\Kelton\analyses\KW006\KW006_07302024_rec_D1_Sub\KW006_07302024_session')
superR(1).root = root;
superS(1).sess = sess;

load('D:\Data\Kelton\analyses\KW006\KW006_07312024_rec_D2_Sub\KW006_07312024_rec_D2_Sub_root')
load('D:\Data\Kelton\analyses\KW006\KW006_07312024_rec_D2_Sub\KW006_07312024_session')
superR(2).root = root;
superS(2).sess = sess;

load('D:\Data\Kelton\analyses\KW007\KW007_08262024_rec_D1_Sub\KW007_08262024_rec_D1_Sub_root')
load('D:\Data\Kelton\analyses\KW007\KW007_08262024_rec_D1_Sub\KW007_08262024_session')
superR(3).root = root;
superS(3).sess = sess;

load('D:\Data\Kelton\analyses\KW007\KW007_08272024_rec_D2_Sub\KW007_08272024_rec_D2_Sub_root')
load('D:\Data\Kelton\analyses\KW007\KW007_08272024_rec_D2_Sub\KW007_08272024_session')
superR(4).root = root;
superS(4).sess = sess;

load('D:\Data\Kelton\analyses\KW008\KW008_08142024_rec_D1_Sub\KW008_08142024_rec_D1_Sub_root')
load('D:\Data\Kelton\analyses\KW008\KW008_08142024_rec_D1_Sub\KW008_08142024_session')
superR(5).root = root;
superS(5).sess = sess;

load('D:\Data\Kelton\analyses\KW008\KW008_08152024_rec_D2_Sub\KW008_08152024_rec_D2_Sub_root')
load('D:\Data\Kelton\analyses\KW008\KW008_08152024_rec_D2_Sub\KW008_08152024_session')
superR(6).root = root;
superS(6).sess = sess;

%% Summarize good/mua/noise per region per session

cd(sdir)

nSess = length(superR);
unitNames = ["Good","MUA","Noise"];
regionNames = ["Cortex","Subiculum","Dentate Gyrus"];

for i = 1:nSess
    unitbreakdown(i,:) = [numel(superR(i).root.good) numel(superR(i).root.mua) numel(superR(i).root.noise)];
    regionbreakdown(i,:) = [numel(superR(i).root.ccs_g_ctx) numel(superR(i).root.ccs_g_sub) numel(superR(i).root.ccs_g_dg)];
    
    pieUnits = figure;
    piechart(unitbreakdown(i,:),unitNames,"LabelStyle","namedata");
    title(superR(i).root.name)
    sname = [superR(i).root.name '_unitBreakdown'];
    saveas(pieUnits,sname,'png')
    
    pieRegions = figure;
    piechart(regionbreakdown(i,:),regionNames,"LabelStyle","namedata");
    colororder reef
    title(["Good Units ", superR(i).root.name])
    sname = [superR(i).root.name '_regionBreakdownGood'];
    saveas(pieRegions,sname,'png')

    close all
end
%% Overall averages across sessions, mice
unitMean = mean(unitbreakdown);
unitStd  = std(unitbreakdown);
unitLoc  = [1 2 3];

figure; hold on;  
b1 = bar(mean(unitbreakdown),'facecolor', 'flat');
errorbar(unitLoc,unitMean,unitStd,'k.')
b1.CData(1,:) = [0 0.4470 0.7410];
b1.CData(2,:) = [0.8500 0.3250 0.0980];
b1.CData(3,:) = [0.9290 0.6940 0.1250];
xs = rand(nSess,1) * 0.2;
scatter(1+xs,unitbreakdown(:,1),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
scatter(2+xs,unitbreakdown(:,2),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
scatter(3+xs,unitbreakdown(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[1 0.75 0.25])
xticks(unitLoc)
xticklabels(unitNames)
ylabel('Unit Count')
set(gca,'FontSize',12,'FontName','Arial')

regionMean = mean(regionbreakdown);
regionStd  = std(regionbreakdown);
regionLoc  = [1 2 3];

figure; hold on;  
b2 = bar(mean(regionbreakdown),'facecolor', 'flat');
errorbar(regionLoc,regionMean,regionStd,'k.')
b2.CData(1,:) = [0.8660    0.3290         0];
b2.CData(2,:) = [0.3290    0.7130    1.0000];
b2.CData(3,:) = [0.0660    0.4430    0.7450];
xs = rand(nSess,1) * 0.2;
scatter(1+xs,regionbreakdown(:,1),'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
scatter(2+xs,regionbreakdown(:,2),'MarkerEdgeColor','k','MarkerFaceColor',[0.25 0.25 0.75])
scatter(3+xs,regionbreakdown(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.25])
xticks(regionLoc)
xticklabels(regionNames)
ylabel('Unit Count')
set(gca,'FontSize',12,'FontName','Arial')

%% Subiculum Stats per session

for i = 1:nSess
    subInds = superR(i).root.ccinds_g_sub;
    subCCs  = superR(i).root.ccs_g_sub;
    nSub = length(subInds);
    for j = 1:nSub
        presUnits(i).present(j) = presenceDelta(superR(i).root,superS(i).sess,subCCs(j));
        [tmpedges,tmpfr,tmpmdl] = plot_frXvel(superR(i).root,subCCs(j),superS(i).sess,2,0);
        velUnitsP(i).velocity(j) = tmpmdl.p < 0.05;
        velUnitsR(i).slope(j) = tmpmdl.b;
        frAll(i).fr(j) = superR(i).root.info.fr(subInds(j));
        frRun(i).fr(j) = mean(tmpfr);
    end
end

%% Summarize stats in Subiculum
velPos = [];
frAllCat = [];
frRunCat = [];

for i = 1:nSess
    presProb(i) = sum(presUnits(i).present)/numel(presUnits(i).present);
    tmppres = logical(presUnits(i).present);
    nUnits(i) = sum(tmppres);
    velProb(i) = mean(velUnitsP(i).velocity(tmppres));
    sigUnits = velUnitsP(i).velocity(tmppres);   %Only significant units
    velSlope = velUnitsR(i).slope(tmppres) > 0; %Positive Speed mod
    velPos = [velPos; mean(velSlope(sigUnits))];    % Percentage +Vel
    velNeg(i) = 1 - velPos(i);  % Percentage -Vel
    frAllCat = [frAllCat; frAll(i).fr(tmppres)'];
    frRunCat = [frRunCat; frRun(i).fr(tmppres)'];
end

%% Plot FR Stats in Subiculum

figure; hold on
for i = 1:nSess
    tmppres = logical(presUnits(i).present);
    errorbar(i,mean(frAll(i).fr(tmppres)),std(frAll(i).fr(tmppres)),'k.')
    xs = rand(numel(frAll(i).fr(tmppres)),1) *0.3;
    scatter(i+xs+0.15, frAll(i).fr(tmppres))
end
xlim([0.5 nSess + .5])
xticklabels(["KW006 D1","KW006 D2","KW007 D1","KW007 D2","KW008 D1","KW008 D2"])
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')

figure; hold on
xs1 = rand(numel(frAllCat),1) * 0.3;
boxplot([frAllCat frRunCat])
scatter(1.15+xs1,frAllCat)
scatter(2.15+xs1,frRunCat)
xlim([0.5 2.5])
xticklabels(["Total","Run"])
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')

%% Plot Speed Stats in Subiculum

figure; hold on
b1 = bar(velProb,'BarWidth',0.5);
ylabel('Velocity Tuning (Probability)')
ylim([0 1])
xlim([0.5 nSess + .5])
xticklabels(["KW006 D1","KW006 D2","KW007 D1","KW007 D2","KW008 D1","KW008 D2"])
set(gca,'FontSize',12,'FontName','Arial')

figure; hold on
bar(mean(velProb),'BarWidth',0.5)
xs = rand(nSess,1) * 0.2;
scatter(1+xs,velProb,'filled')
ylim([0 1])
ylabel('Velocity Tuning (Probability)')
xlim([0.5 1.5])
xticks([])
set(gca,'FontSize',12,'FontName','Arial')

figure; hold on
b1 = bar([velPos, velNeg'],'stacked','BarWidth',0.5);
ylabel('Velocity Tuning (Probability)')
legend({'Positive','Negative'})
xlim([0.5 nSess + .5])
xticklabels(["KW006 D1","KW006 D2","KW007 D1","KW007 D2","KW008 D1","KW008 D2"])
set(gca,'FontSize',12,'FontName','Arial')

figure; hold on
bar([mean(velPos) mean(velNeg)],'BarWidth',0.5)
xs = rand(nSess,1) * 0.2;
scatter(1+xs,velPos,'filled')
scatter(2+xs,velNeg,'filled')
ylim([0 1])
ylabel('Velocity Tuning (Probability)')
xlim([0.5 2.5])
xticklabels({"Positive","Negative"})
set(gca,'FontSize',12,'FontName','Arial')

for i = 1:nSess
    tmppres = logical(presUnits(i).present);
    errorbar(i,mean(frAll(i).fr(tmppres)),std(frAll(i).fr(tmppres)),'k.')
    xs = rand(numel(frAll(i).fr(tmppres)),1) *0.3;
    scatter(i+xs+0.15, frAll(i).fr(tmppres))
end
xlim([0.5 nSess + .5])
xticklabels(["KW006 D1","KW006 D2","KW007 D1","KW007 D2","KW008 D1","KW008 D2"])
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')

