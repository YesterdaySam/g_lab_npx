% Summarize behavior script

mousedir = 'D:\Data\Kelton\analyses\KW006';

overwriteFlag = 0;

cd(mousedir)

dirlist = dir;
dirlist = dirlist([dirlist.isdir]);

% Pre Gen Data Table
varTypes = ["string","string","double","double","double","double","double"];
varNames = ["SessName","SessType","Day","Length(min)","nLaps","uLicks","uVelocity"];
sz = [1 length(varTypes)];
sumT = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

rowCt = 1;

for i = 1:length(dirlist)
    cd(fullfile(dirlist(i).folder,dirlist(i).name));

    try 
        tmpSess = dir("*session.mat");
        load(tmpSess.name)
    catch
        cd(mousedir)
        continue
    end
    sumT.SessName(rowCt)        = sess.name;
    if ~isempty(strfind(sess.name,'training'))
        sumT.SessType(rowCt)    = 'training';
    else
        sumT.SessType(rowCt)    = 'recording';
    end
    tmpStr                      = strsplit(sess.name,'_D');
    sumT.Day(rowCt)             = str2double(strtok(tmpStr{2},'_'));
    sumT.("Length(min)")(rowCt) = sess.ts(end)/60;
    sumT.nLaps(rowCt)           = sess.nlaps;
    sumT.uLicks(rowCt)          = length(sess.lckind) / sess.nlaps; %in licks/lap
    sumT.uVelocity(rowCt)       = mean(sess.velshft);    %in cm/s
    sumT.uVelStp(rowCt)         = mean(sess.velshft(sess.velshft > 0 & sess.velshft < 2)); %in cm/s, low run speed velocity
    sumT.uVelRun(rowCt)         = mean(sess.velshft(sess.velshft > 2)); %in cm/s, high run speed velocity
    sumT.totDay(rowCt)          = rowCt;
    rowCt = rowCt + 1;
end

cd(mousedir)
save('bhvrSummaryT','sumT')

disp(['Finished behavior summary for ' mousedir(end-4:end)])

%% Summary graphics

cd(mousedir)

% Training plots
trnInds = sumT.SessType == 'training';
recInds = sumT.SessType == 'recording';

figure;
set(gcf,'units','normalized','position',[0.4 0.15 0.5 0.7])

subplot(2,2,1); hold on
plot(sumT.Day(trnInds), sumT.("Length(min)")(trnInds),'k','LineWidth',2)
plot(sumT.totDay(recInds), sumT.("Length(min)")(recInds),'b','LineWidth',2)
legend({'training','recording'},'location','se')
ylabel('Recording Length (min)')
ylim([0 inf])

subplot(2,2,2); hold on
plot(sumT.Day(trnInds), sumT.nLaps(trnInds),'k','LineWidth',2)
plot(sumT.totDay(recInds), sumT.nLaps(recInds),'b','LineWidth',2)
ylabel('Number of Laps')
ylim([0 inf])

subplot(2,2,3); hold on
plot(sumT.Day(trnInds),sumT.uLicks(trnInds),'k','LineWidth',2)
plot(sumT.totDay(recInds),sumT.uLicks(recInds),'b','LineWidth',2)
xlabel('Behavior Day')
ylabel('Average Licks/lap')
ylim([0 inf])

subplot(2,2,4); hold on
% plot(sumT.Day(trnInds),sumT.uVelStp(trnInds),'k','LineWidth',2)
plot(sumT.Day(trnInds),sumT.uVelRun(trnInds),'k','LineWidth',2)
% plot(sumT.totDay(recInds),sumT.uVelStp(recInds),'b','LineWidth',2)
plot(sumT.totDay(recInds),sumT.uVelRun(recInds),'b','LineWidth',2)
xlabel('Behavior Day')
ylabel('Average Velocity (cm/s)')
ylim([0 inf])

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12);
sgtitle(mousedir(end-4:end))

saveas(gcf,'bhvrSummary1.png')