%% Prep Sessions/Animals/variables for later blocks
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%
% ========================================================================%

datT = import_xldat("D:\Data\Kelton\analyses\group_analyses\Subiculum_DCZ","dat_include_zm.xlsx");
groupSDir = 'D:\Data\Kelton\analyses\group_analyses\Subiculum_DCZ\analysis';
cd(groupSDir) 

mInclude = {'ZM040','ZM042'};

sessType = 2;
for i = 1:height(datT)
    useInds(i) = logical(sum(strcmp(datT.mouse(i),mInclude))) & logical(sum(ismember(sessType,datT.sess_type(i)))); 
end
datT(~useInds,:) = [];  %Clean excluded sessions

saveFlag = 1;
sbase = 'subDCZShift_';
bvName = [sbase 'groupDat_bhvr'];

dbnsz = 0.05;
histoBnsz = 5;
binedges = 0:5:185;
binpos = 0.025:dbnsz:1.825;
wlen = 150;
r1pos = 0.4;    % 10 cm
r2pos = 1.3;      % 100cm
fvncols = [.35 .35 .35; 1 .25 .25]; % Gray vs red
lnlcols = [0.0508 0.4883 0.5273; 0.7852 0.6055 0.2188]; % A teal vs B brown
% rgncols = [0.9961 0.7305 0.4336; 0.3672 0.2969 0.3711]; % Sub gold vs CA1 dull purple

clear ps stats

%% Combine Sub RZ Shift Behavior data
combine_bhvrDat(datT,bvName,groupSDir,2);

%% Load saved behavior data
cd(groupSDir)
load(bvName)
nMice = length(unique(bhvID(:,1)));
mID = unique(bhvID(:,1),'stable');

mA = [42];
mB = [40];

aInd = [];
bInd = [];

for i = 1:size(bhvID,1)
    if ~isempty(find(mA == bhvID(i,1), 1))
        aInd = [aInd; find(bhvID(:,1) == bhvID(i,1))];
    else
        bInd = [bInd; find(bhvID(:,1) == bhvID(i,1))];
    end
end

bvgrp(1).bvInd = aInd;
bvgrp(2).bvInd = bInd;
bvgrp(1).grpname = 'virusA';
bvgrp(2).grpname = 'virusB';
bvgrp(1).mID = bhvID(bvgrp(1).bvInd);
bvgrp(2).mID = bhvID(bvgrp(2).bvInd);
bvgrp(1).n = length(bvgrp(1).bvInd);
bvgrp(2).n = length(bvgrp(2).bvInd);

%% Behavior comparisons
% Learners: P(Nov Rwd) > 0.55

for i = 1:nMice
    uLapRwd50(i,:) = [mean(bvDat(i).preLapRwd(end-49:end)) mean(bvDat(i).pstLapRwd(1:50))];
end
% uLDI = [vertcat(bvDat.uPreLckDI), vertcat(bvDat.uPstLckDI)];
uPsv = [vertcat(bvDat.uPreLckPsv), vertcat(bvDat.uPstLckPsv)];
uLapRwd = [vertcat(bvDat.uPreLapRwd), vertcat(bvDat.uPstLapRwd)];
nLaps = [vertcat(bvDat.preNLap), vertcat(bvDat.pstNLap)];

psvSplitF = plot_2d_bhvr(uPsv,bvgrp(1).bvInd,bvgrp(2).bvInd,lnlcols);
% plot([-1 1],[-0.15 -0.15],'k--')
xlabel('LDI F'); xlim([-1 1])
ylabel('LDI N'); ylim([-1 1])

rwdSplitF = plot_2d_bhvr(uLapRwd,bvgrp(1).bvInd,bvgrp(2).bvInd,lnlcols);
plot([0 1],[0.55 0.55],'k--')
xlabel('P(F Lap Rewarded)'); xlim([0 1])
ylabel('P(N Lap Rewarded)'); ylim([0 1])
legend({'A','B'}, 'Location','sw')

lapSplitF = plot_2d_bhvr(nLaps,bvgrp(1).bvInd,bvgrp(2).bvInd,lnlcols);
xlabel('# Laps F'); xlim([0 140])
ylabel('# Laps N'); ylim([0 140])

if saveFlag
    fsave(psvSplitF,[sbase 'bhv_psvPrePst'],1,1);
    fsave(rwdSplitF,[sbase 'bhv_rwdPrePst'],1,1);
    fsave(lapSplitF,[sbase 'bhv_lapPrePst'],1,1);
end


%% Functions 
function [fhandle] = plot_2d_bhvr(dat,learn,nonlearn,cols)
fhandle = figure; hold on; axis square
set(gcf,'Units','normalized','Position',[0.1 0.4 0.3333 0.1786])
plot(dat(learn,1),dat(learn,2),'o','Color',cols(1,:))
plot(dat(nonlearn,1),dat(nonlearn,2),'o','Color',cols(2,:))
set(gca,'FontSize',16,'FontName','Arial')
end
