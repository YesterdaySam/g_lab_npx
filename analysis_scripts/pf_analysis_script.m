% root = root4;
% sess = sess4;
cd('D:\Data\Kelton\analyses\KW022\KW022_12172024_rec_D2_RMed1')
saveFlag = 0;

metafile = dir('*_meta.mat');
load(metafile.name)

%%
SI = zeros(1,length(root.good));
uFR = zeros(1,length(root.good));
pkFR = zeros(1,length(root.good));

nUnits = length(root.good);

for i = 1:nUnits
    cc = root.good(i);

    % [SI(i),uFR(i),pkFR(i)] = get_SI(root,cc,sess);
    [SI(i),pkFR(i),uFR(i)] = get_PF(root,cc,sess);
end

%% Find good PFs and Plot waterfalls of place fields
sithresh = 0.2;
frthresh = 2;
uThresh  = 10;

hiSI = SI > sithresh;
hiFR = pkFR > frthresh;
lowuFR = uFR < uThresh;

% Find depth and shankID of good fields
goodPFs = root.good(hiSI & hiFR & lowuFR);
clear tmprgn

for i = 1:numel(goodPFs)
    [goodPFs(i,2), goodPFs(i,3), tmprgn(i)] = getDepthByID(root,goodPFs(i));
end
[~,goodPFsSort] = sort(goodPFs(:,3));
goodPFsSort = goodPFs(goodPFsSort,:);

tmpfig = plot_unitsXpos(root,sess,goodPFs);
if meta.prbType == 'NPX2.0'
    tmpfig2 = figure; histogram(goodPFs(:,3))
    set(gcf,'units','normalized','position',[0.2 0.4 0.2 0.4])
    xlabel('Shank ID')
    ylabel('"Good" PF count')
    set(gca,'FontSize',12,'FontName','Arial')
end

if saveFlag
    saveas(tmpfig,'goodSIFR_waterfall','png');
    saveas(tmpfig2,'goodSIFRxShank','png');
    save([root.name '_goodPFs'],'goodPFs')
end

%% Find good PFs by region

goodCA1 = goodPFs(tmprgn == "CA1",:);
goodSub = goodPFs(tmprgn == "Sub",:);

ca1Waterfall = plot_unitsXpos(root,sess,goodCA1);
subWaterfall = plot_unitsXpos(root,sess,goodSub);

allCA1 = sum(root.info.region(root.goodind) == "CA1");
allSub = sum(root.info.region(root.goodind) == "Sub");

prGoodPFsub = length(goodSub)/allSub;
prGoodPFca1 = length(goodCA1)/allCA1;

regionGoodPFpr = figure;
set(gcf,'units','normalized','position',[0.45 0.3 0.1 0.4])
bar([prGoodPFsub prGoodPFca1])
xticklabels({"Sub" "CA1"})
ylabel('"Good" PF Probability')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(ca1Waterfall,'goodCA1_waterfall','png');
    saveas(subWaterfall,'goodSub_waterfall','png');
    saveas(regionGoodPFpr,'goodPF_regionPr','png');
end

%% Find speed coding units

ps = zeros(1,length(root.good));
rs = zeros(1,length(root.good));
bs = zeros(1,length(root.good));

nUnits = length(root.good);

for i = 1:nUnits
    cc = root.good(i);

    [~,~,tmpModel] = plot_frXvel(root,cc,sess,2,0);

    ps(i) = tmpModel.p;
    rs(i) = tmpModel.r;
    bs(i) = tmpModel.b;
end

psthresh = 0.05;
rsthresh = 0.4;

lops = ps < psthresh;
hirs = rs > rsthresh;

% Find depth and shankID of good fields
goodVels = root.good(lops & hirs);
goodbs   = bs(lops & hirs)';

for i = 1:numel(goodVels)
    [goodVels(i,2), goodVels(i,3), tmprgn(i)] = getDepthByID(root,goodVels(i));
end

if meta.prbType == 'NPX2.0'
    tmpfigV = figure; histogram(goodVels(:,3))
    set(gcf,'units','normalized','position',[0.2 0.4 0.2 0.4])
    xlabel('Shank ID')
    ylabel('"Good" Velocity count')
    set(gca,'FontSize',12,'FontName','Arial')
end

if saveFlag
    saveas(tmpfigV,'goodVelxShank','png');
    save([root.name '_goodVels'],'goodVels')
end

%% Find good speed units by region

goodCA1 = goodVels(tmprgn == "CA1",:);
goodSub = goodVels(tmprgn == "Sub",:);
goodCA1bs = goodbs(tmprgn == "CA1",:);
goodSubbs = goodbs(tmprgn == "Sub",:);

allCA1 = sum(root.info.region(root.goodind) == "CA1");
allSub = sum(root.info.region(root.goodind) == "Sub");

prGoodPFsub = length(goodSub)/allSub;
prGoodPFca1 = length(goodCA1)/allCA1;

regionGoodVpr = figure;
set(gcf,'units','normalized','position',[0.45 0.3 0.1 0.4])
bar([prGoodPFsub prGoodPFca1])
xticklabels({"Sub" "CA1"})
ylim([0 1])
ylabel('"Good" Speed Cell Prob.')
set(gca,'FontSize',12,'FontName','Arial')

ca1Slopes = [sum(goodCA1bs > 0)/length(goodCA1) sum(goodCA1bs < 0)/length(goodCA1)];
subSlopes = [sum(goodSubbs > 0)/length(goodSub) sum(goodSubbs < 0)/length(goodSub)];

regionSlopeProb = figure;
set(gcf,'units','normalized','position',[0.45 0.3 0.1 0.4])
bar([subSlopes; ca1Slopes])
xticklabels({"Sub" "CA1"})
ylim([0 1])
legend({'Positve','Negative'},'Location','northeast')
ylabel('Probability')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(regionGoodVpr,'goodVel_regionPr','png');
    saveas(regionSlopeProb,'goodVel_slopePr','png');
end

%% compare shifted spike trains for pfs
cc = 321;
nShufs = 1000;

[tmpSI,tmppkFR,tmpuFR,~,tmpbinFR,tmpedges] = get_PF(root,cc,sess);
shufSI = zeros(nShufs,1);
shufuFR = zeros(nShufs,1);
shufpkFR = zeros(nShufs,1);
shufbinFR = zeros(nShufs,length(tmpedges)-1);

for i = 1:nShufs
    shiftroot = shiftTrain(root,sess);
    % [shufSI(i),shufpkFR(i),shufuFR(i),~,shufbinFR(i,:)] = get_PF(shiftroot,cc,sess); 
    [shufSI(i),shufuFR(i),shufpkFR(i),~,~,shufbinFR(i,:)] = get_SI(shiftroot,cc,sess);
end

sem = rmmissing(std(shufbinFR,'omitnan')/sqrt(nShufs));
ciup = mean(shufbinFR,'omitnan') + sem*1.96;
cidn = mean(shufbinFR,'omitnan') - sem*1.96;

figure; hold on
plot(tmpedges(1:end-1)*100,tmpbinFR,'b')
patch(100*[tmpedges(1:length(tmpbinFR)),fliplr(tmpedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
plot(tmpedges(1:end-1)*100,mean(shufbinFR),'k')
xlabel('Position (cm)')
ylabel('Firing Rate')
title(['Cell ' num2str(cc)])
legend({"True","Avg Shuffle"})
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

figure; hold on
histogram(shufSI,30)
if tmpSI > prctile(shufSI,99)
    plot([tmpSI, tmpSI], [0 100], 'r--')
else
    plot([tmpSI, tmpSI], [0 100], 'k--')
end
xlabel('Spatial Info.')
ylabel('Count')
title(['Cell ' num2str(cc)])
legend({"Shuffle","True"})
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

shufTest = tmpSI > prctile(shufSI,99)

%% Batch shuffle all good units
nShufs = 1000;

shufTest = zeros(numel(root.good),1);

for i = 1:numel(root.good)
    cc = root.good(i);

    tmpSI = get_SI(root,cc,sess);
    shufSI = zeros(nShufs,1);

    for j = 1:nShufs
        shiftroot = shiftTrain(root,sess);
        [shufSI(j)] = get_SI(shiftroot,cc,sess);
    end

    shufTest(i) = tmpSI > prctile(shufSI,99);
end

figure; imagesc(shufTest)