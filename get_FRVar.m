function [root] = get_FRVar(root,sess,plotflag)

arguments
    root
    sess
    plotflag = 1
end

sdMu = [];
for i = 1:height(root.info)
    unit = root.info.cluster_id(i);
    [~,~,tmpSDMu] = get_presence(root,unit,sess);
    sdMu = [sdMu; tmpSDMu];
end

root.info.frVar = sdMu;

if plotflag
    figure; hold on
    plot(root.info.cluster_id(root.goodind),root.info.frVar(root.goodind),'b')
    plot(root.info.cluster_id(root.muaind),root.info.frVar(root.muaind),'r')
    plot(root.info.cluster_id(root.noiseind),root.info.frVar(root.noiseind),'Color',[.6 .6 .6])
    plot([1 max(root.info.cluster_id)], [1.4 1.4],'k--')
    xlabel('Cluster ID')
    ylabel('FR Stdev/Mean')

    hiVars = find(root.goodind & root.info.frVar > 1.4);
    nSubFigs = ceil(length(hiVars)/16);
    ct = 0;
    for i = 1:nSubFigs
        figure
        for j = 1:16
            try
                subplot(4,4,j)
                hold on
                unit = find(root.info.cluster_id == hiVars(j+ct));
                [tmpCts] = get_presence(root,unit,sess);
                bar(tmpCts)
                plot([0 length(tmpCts)],[mean(tmpCts) mean(tmpCts)], 'k--')
                xlim([0 length(tmpCts)])
                title(['Unit ' num2str(root.info.cluster_id(hiVars(j+ct)))])
                legend({['FR Var' num2str(root.info.frVar(hiVars(j+ct)))]},'Location','se','FontSize',8)
            catch
                continue
            end
        end
        ct = ct + 16;
    end
end

end