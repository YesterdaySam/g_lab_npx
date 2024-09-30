% Summarize behavior script across many mice
% Run import_script.m and bhvrSummary.m scripts first to generate behavior
% tables prior to running this

micedir = 'D:\Data\Kelton\analyses';

cd(micedir)

dirlist = dir;
dirlist = dirlist([dirlist.isdir]);

superT = struct;
mNames = {};
rowCt = 1;

for i = 1:length(dirlist)
    cd(fullfile(dirlist(i).folder,dirlist(i).name));

    try 
        tmpT = dir("*SummaryT.mat");
        load(tmpT.name)
    catch
        cd(micedir)
        continue
    end
    mNames{i} = sumT.SessName{1}(1:5);
    sumT.mouse = repmat(mNames{i},height(sumT),1);
    sumT.LapsPerMin = sumT.nLaps ./ sumT.("Length(min)");
    superT(rowCt).sumT = sumT;
    rowCt = rowCt + 1;

end

% superT.LapsPerMin = superT.nLaps ./ superT.("Length(min)");

%% Prep for plotting
sdir = 'D:\Data\Kelton\analyses\group_analyses\KW001-12_behavior';

cd(sdir)

nMice = length(superT);
% cmap1 = winter(nMice);
cmap2 = cool(nMice);

%% Rec Length

figure;
set(gcf,'units','normalized','position',[0.3 0.15 0.25 0.6])

for i = 1:nMice
    sumT = superT(i).sumT;

    trnInds = sumT.SessType == 'training';
    recInds = sumT.SessType == 'recording';

    subplot(2,1,1); hold on
    plot(sumT.Day(trnInds), sumT.("Length(min)")(trnInds),'Color',cmap2(i,:),'LineWidth',1)
    legCell1{i} = sumT.mouse(1,:);
    ylabel('Recording Length (min)')
    % ylim([0 inf])
    ys = ylim;

    subplot(2,1,2); hold on
    plot(sumT.Day(recInds), sumT.("Length(min)")(recInds),'Color',cmap2(i,:),'LineWidth',1)
    ylabel('Recording Length (min)')
    ylim(ys)
    xlim([0 inf])
end
subplot(2,1,1); hold on
title('Training')
% legend(legCell1,'location','se')
% xs = xlim;
% xlim([xs(1) xs(2) + 10])

subplot(2,1,2); hold on
title('Experiment')
xlim([0 10])
xlabel('Behavior Day')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12);

saveas(gcf,'recLenSummary.png')

%% nLaps

figure;
set(gcf,'units','normalized','position',[0.3 0.15 0.25 0.6])

for i = 1:nMice
    sumT = superT(i).sumT;

    trnInds = sumT.SessType == 'training';
    recInds = sumT.SessType == 'recording';

    subplot(2,1,1); hold on
    plot(sumT.Day(trnInds), sumT.nLaps(trnInds),'Color',cmap2(i,:),'LineWidth',1)
    legCell1{i} = sumT.mouse(1,:);
    ylabel('Number of Laps')
    % ylim([0 inf])
    ys = ylim;

    subplot(2,1,2); hold on
    plot(sumT.Day(recInds), sumT.nLaps(recInds),'Color',cmap2(i,:),'LineWidth',1)
    ylabel('Number of Laps')
    ylim(ys)
    xlim([0 inf])
end
subplot(2,1,1); hold on
title('Training')

subplot(2,1,2); hold on
title('Experiment')
xlim([0 10])
xlabel('Behavior Day')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12);

saveas(gcf,'nLapsSummary.png')

%% Average licks/lap

figure;
set(gcf,'units','normalized','position',[0.3 0.15 0.25 0.6])

for i = 1:nMice
    sumT = superT(i).sumT;

    trnInds = sumT.SessType == 'training';
    recInds = sumT.SessType == 'recording';

    subplot(2,1,1); hold on
    plot(sumT.Day(trnInds), sumT.uLicks(trnInds),'Color',cmap2(i,:),'LineWidth',1)
    legCell1{i} = sumT.mouse(1,:);
    ylabel('Average Licks/lap')
    ylim([0 400])
    ys = ylim;

    subplot(2,1,2); hold on
    plot(sumT.Day(recInds), sumT.uLicks(recInds),'Color',cmap2(i,:),'LineWidth',1)
    ylabel('Average Licks/lap')
    % ylim(ys)
    xlim([0 inf])
end
subplot(2,1,1); hold on
title('Training')
legend(legCell1,'location','ne')
% xs = xlim;
% xlim([xs(1) xs(2) + 10])

subplot(2,1,2); hold on
title('Experiment')
xlim([0 10])
xlabel('Behavior Day')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12);

saveas(gcf,'uLicksSummary.png')

%% Average run velocity

figure;
set(gcf,'units','normalized','position',[0.3 0.15 0.25 0.6])

for i = 1:nMice
    sumT = superT(i).sumT;

    trnInds = sumT.SessType == 'training';
    recInds = sumT.SessType == 'recording';

    subplot(2,1,1); hold on
    plot(sumT.Day(trnInds), sumT.uVelRun(trnInds),'Color',cmap2(i,:),'LineWidth',1)
    legCell1{i} = sumT.mouse(1,:);
    ylabel('Average Velocity (cm/s)')
    % ylim([0 400])
    ys = ylim;

    subplot(2,1,2); hold on
    plot(sumT.Day(recInds), sumT.uVelRun(recInds),'Color',cmap2(i,:),'LineWidth',1)
    ylabel('Average Velocity (cm/s)')
    ylim(ys)
    xlim([0 inf])
end
subplot(2,1,1); hold on
title('Training')
% legend(legCell1,'location','ne')
% xs = xlim;
% xlim([xs(1) xs(2) + 10])

subplot(2,1,2); hold on
title('Experiment')
xlim([0 10])
xlabel('Behavior Day')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12);

saveas(gcf,'uVelRunSummary.png')

%% Average Laps/min

figure;
set(gcf,'units','normalized','position',[0.3 0.15 0.25 0.6])

for i = 1:nMice
    sumT = superT(i).sumT;

    trnInds = sumT.SessType == 'training';
    recInds = sumT.SessType == 'recording';

    subplot(2,1,1); hold on
    plot(sumT.Day(trnInds), sumT.LapsPerMin(trnInds),'Color',cmap2(i,:),'LineWidth',1)
    legCell1{i} = sumT.mouse(1,:);
    ylabel('Laps/min')
    % ylim([0 400])
    ys = ylim;

    subplot(2,1,2); hold on
    plot(sumT.Day(recInds), sumT.LapsPerMin(recInds),'Color',cmap2(i,:),'LineWidth',1)
    ylabel('Laps/min')
    ylim(ys)
    xlim([0 inf])
end
subplot(2,1,1); hold on
title('Training')
% legend(legCell1,'location','ne')
% xs = xlim;
% xlim([xs(1) xs(2) + 10])

subplot(2,1,2); hold on
title('Experiment')
xlim([0 10])
xlabel('Behavior Day')

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12);

saveas(gcf,'lapsPerMinSummary.png')
