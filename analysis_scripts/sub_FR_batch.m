%% Prepare data table and run parameters
parentDir = "D:\Data\Kelton\analyses\group_analyses\Subiculum_FR"; 
subDatT = import_xldat(parentDir,"dat_include.xlsx");
cd(parentDir) 

cleanInds = subDatT.include == 0;
subDatT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;

%% Load in data

frAll = [];
frRun = [];
frSit = [];
mouse = [];
recID = [];

for i = 1:height(subDatT)
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
        root = get_estCellType(root,15,5,100,1);    %Screen putative interneurons based on width and FWHM thresholds, not FR
    end
    root = get_FRVar(root,sess,0);  %Screen for overly high variance of firing rate

    for j = 1:nShank
        disp(['Shank' num2str(j-1)])
        if subDatT{i,6+j}{1} == 'sub'
            lyrShUnits = find(root.info.shankID == j-1 & root.info.lyrID == 1 & root.goodind & root.info.uType == 1 & root.info.frVar < 1.4);
            disp(['N = ' num2str(sum(root.info.uType(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1))),...
                ' good pyrs and N = ' num2str(sum(~root.info.uType(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1))) ' good INs'])
            disp(['N = ' num2str(sum(root.info.frVar(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1) < 1.4)),...
                ' stable units and N = ' num2str(sum(root.info.frVar(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1)>=1.4)) ' high variance units'])
            frAll = [frAll; root.info.fr(lyrShUnits)];
            for k = 1:length(lyrShUnits)
                cc = root.info.cluster_id(lyrShUnits(k));
                tmpFR = get_frStandVRun(root,cc,sess,4);
                frSit = [frSit; tmpFR(:,1)];
                frRun = [frRun; tmpFR(:,2)];
            end
            mouse = [mouse; repmat(str2num(subDatT.mouse{i}(end-1:end)), length(lyrShUnits),1)];
            recID = [recID; ones(length(lyrShUnits),1)*i];
        end
    end
end

cd(parentDir)

%% Observational statistics

mNames = unique(mouse);
nMice = numel(unique(mouse));
nSess = max(recID);
nUnits = length(frSit);

for i = 1:nMice
    mouseInds = mouse == mNames(i);
    unitsByMouse(i) = sum(mouseInds);
    recsByMouse(i) = numel(unique(recID(mouseInds)));
end
for i = 1:nSess
    recInds = recID == i;
    unitsByRec(i) = sum(recInds);
end

disp(['Total units: ', num2str(nUnits)])
disp(['Number of mice: ', num2str(nMice)])
disp(['Number of recordings: ', num2str(nSess)])
disp(['Number of recordings per mouse: ', num2str(mean(recsByMouse))])
disp(['Mean units per mouse: ', num2str(mean(unitsByMouse)), '; Stdev units per mouse: ', num2str(std(unitsByMouse))])
disp(['Mean units per recording: ', num2str(mean(unitsByRec)), '; Stdev units per recording: ', num2str(std(unitsByRec))])
disp(['Mean FR standing: ', num2str(mean(frSit)), '; Stdev FR standing: ', num2str(std(frSit))])
disp(['Mean FR running: ', num2str(mean(frRun)), '; Stdev FR running: ', num2str(std(frRun))])
disp(['Mean FR overall: ', num2str(mean(frAll)), '; Stdev FR overall: ', num2str(std(frAll))])

%% Statistical tests

ps = [];
stats = [];

[h, ps.kstest, stats.kstetst] = kstest([frSit; frRun]);

if h
    [ps.wilcRM,~,stats.wilcRM] = signrank(frSit,frRun);
else
    [~,ps.ttestRM,~,stats.ttestRM] = ttest2(frSit, frRun);
end

%% Violinplot

catDat = [frAll; frSit; frRun];
catLbs = [repmat("Session",nUnits,1); repmat("Standing",nUnits,1); repmat("Running",nUnits,1);];
catOrder = {'Session','Standing','Running'};

violinFig = figure;
violinplot(catDat,catLbs,'GroupOrder',catOrder);
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag; saveas(violinFig,['subFR_FR_distro_violin_stablePyrs.png']); end

%% Boxplot
randXs = -0.025 + 0.05*rand(nUnits,1);

meanFig = figure; hold on
set(gcf,"Position",[680 458 330 420])
boxplot([frSit, frRun], {'Standing','Running'},'Symbol','','Colors','k','BoxStyle','filled','MedianStyle','target')
plot([1.1*ones(nUnits,1)+randXs, 1.9*ones(nUnits,1)+randXs]',[frSit, frRun]','-o','Color',[0.6 0.6 0.6])
plot([1.1 1.9]',[mean(frSit),mean(frRun)]','k-o','LineWidth',3)
xlim([.5 2.5])
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')

if saveFlag; saveas(meanFig,['subFR_FR_distro_means_stablePyrs.png']); end

%% Boxplot
randXs = -0.2 + 0.4*rand(nUnits,1);

boxFig = figure; hold on
set(gcf,"Position",[680 458 330 420])
plot([1*ones(nUnits,1)+randXs, 2*ones(nUnits,1)+randXs, 3*ones(nUnits,1)+randXs]',[frSit, frRun, frAll]','k.') %'Color',[0.6 0.6 0.6])
boxplot([frSit, frRun, frAll], {'Standing','Running','Overall'},'Symbol','','Colors','k')
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')
if saveFlag; saveas(boxFig,['subFR_FR_distro_boxplot_stablePyrs.png']); end

%% Running only Boxplot
randXs = -0.15 + 0.3*rand(nUnits,1);

boxFig2 = figure; hold on
set(gcf,"Position",[680 458 330 420])
plot([1*ones(nUnits,1)+randXs]',[frRun]','k.') %'Color',[0.6 0.6 0.6])
boxplot([frRun], {'Running'},'Symbol','','Colors','k')
ylabel('Firing Rate (Hz)')
set(gca,'FontSize',12,'FontName','Arial')
if saveFlag; saveas(boxFig2,['subFR_FR_distro_boxplot_stablePyrs_runonly.png']); end

%% SPW-R analyses

% for i = 1:height(subDatT)
%     cd(subDatT.fpath{i})
% 
%     rootfile = dir("*_root.mat");
%     load(rootfile.name)
%     sessfile = dir("*_session.mat");
%     load(sessfile.name)
%     disp(root.name)
% 
%     nShank = numel(unique(root.info.shankID));
% 
%     % root = get_layerUnits(root,100);
%     try
%         root.info.uType(1);
%     catch
%         root = get_estCellType(root,15,5,100,1);    %Screen putative interneurons based on width and FWHM thresholds, not FR
%     end
%     root = get_FRVar(root,sess,0);  %Screen for overly high variance of firing rate
% 
%     for j = 1:nShank
%         disp(['Shank' num2str(j-1)])
%         if subDatT{i,6+j}{1} == 'sub'
%             lyrShUnits = find(root.info.shankID == j-1 & root.info.lyrID == 1 & root.goodind & root.info.uType == 1 & root.info.frVar < 1.4);
%             disp(['N = ' num2str(sum(root.info.uType(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1))),...
%                 ' good pyrs and N = ' num2str(sum(~root.info.uType(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1))) ' good INs'])
%             disp(['N = ' num2str(sum(root.info.frVar(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1) < 1.4)),...
%                 ' stable units and N = ' num2str(sum(root.info.frVar(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1)>=1.4)) ' high variance units'])
%             frAll = [frAll; root.info.fr(lyrShUnits)];
%             for k = 1:length(lyrShUnits)
%                 cc = root.info.cluster_id(lyrShUnits(k));
%                 tmpFR = get_frStandVRun(root,cc,sess,4);
%                 frSit = [frSit; tmpFR(:,1)];
%                 frRun = [frRun; tmpFR(:,2)];
%             end
%             mouse = [mouse; repmat(str2num(subDatT.mouse{i}(end-1:end)), length(lyrShUnits),1)];
%             recID = [recID; ones(length(lyrShUnits),1)*i];
%         end
%     end
% end
% 
% cd(parentDir)

