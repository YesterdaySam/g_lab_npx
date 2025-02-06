function [fhandle] = plot_spikeDensity(root,binsize)
%% Plots the binned spike density regardless of cluster across time
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% binsize = size of temporal bins, default = 60 (sec)
%
% Outputs:
% fhandle = handle to figure
%
% Created 2/3/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

nShanks = numel(unique(root.info.shankID));

binedges = (0:root.fspulse*binsize:length(root.syncpulse))/root.fspulse;
tmpMuCount = [];
tmpGdCount = [];
tmpBdCount = [];

for i = 1:max(root.info.shankID)+1
    muIncl = zeros(1,length(root.ts));
    gdIncl = zeros(1,length(root.ts));
    bdIncl = zeros(1,length(root.ts));

    for j = 1:height(root.info)
        if root.info.shankID(j) == i-1 & root.info.group{j}(1) ~= 'm' % If right shank and (m)ua
            muIncl(root.cl == j) = 1;
        end
        if root.info.shankID(j) == i-1 & root.info.group{j}(1) == 'g' % If right shank and (g)ood
            gdIncl(root.cl == j) = 1;
        end
        if root.info.shankID(j) == i-1 & root.info.group{j}(1) == 'n' % If right shank and (n)oise
            bdIncl(root.cl == j) = 1;
        end
    end
    
    tmpMuCount(i,:) = histcounts(root.ts(logical(muIncl)),binedges);
    tmpGdCount(i,:) = histcounts(root.ts(logical(gdIncl)),binedges);
    tmpBdCount(i,:) = histcounts(root.ts(logical(bdIncl)),binedges);

end

%% Plot

fhandle = figure;
set(gcf,'units','normalized','position',[0.4 0.05 0.3 0.85])

cmapGrey = gray(nShanks + 1);
cmapRed  = hot(nShanks + 5);
cmapSky  = sky(nShanks + 2);

subplot(3,1,1); hold on
plot(binedges(2:end),mean(tmpGdCount),'b','LineWidth',2)
for i = 1:max(root.info.shankID)+1
    plot(binedges(2:end),tmpGdCount(i,:),'Color',cmapSky(i,:))
end
title('Good spikes')
set(gca,'FontSize',12,'FontName','Arial')

subplot(3,1,2); hold on
plot(binedges(2:end),mean(tmpMuCount),'r','LineWidth',2)
for i = 1:max(root.info.shankID)+1
    plot(binedges(2:end),tmpMuCount(i,:),'Color',cmapRed(i+1,:))
end
title('MUA spikes')
ylabel('Spike density (counts)')
set(gca,'FontSize',12,'FontName','Arial')

subplot(3,1,3); hold on
plot(binedges(2:end),mean(tmpBdCount),'k','LineWidth',2)
for i = 1:max(root.info.shankID)+1
    plot(binedges(2:end),tmpBdCount(i,:),'Color',cmapGrey(i,:))
end
title('Noise spikes')
xlabel('Time (sec)')
set(gca,'FontSize',12,'FontName','Arial')

end