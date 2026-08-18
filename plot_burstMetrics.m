function [fhandle] = plot_burstMetrics(root,sess,useCC,rgdat)

for i = 1:length(root.good)
    cc = root.good(i);
    [burstIndex(i),burstISI(i),burstLens(i)] = get_burstIndex(root,sess,cc);
end

nrgns = unique(rgdat);
% linespec = {'co','mo','yo','ko'};
fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
fhandle = fixRatio(fhandle);
for i = 1:length(nrgns)
    rgn = rgdat(root.goodind) == nrgns(i);
    if nrgns(i) == 1    % CA1
        % linespec = [0.2891 0.1367 0.4648]; % Dark purple
        linespec = [0.7852 0.6055 0.2188]; % Med brown
    elseif nrgns(i) == 2    % Sub
        linespec = [0.0508 0.4883 0.5273]; % Dark teal
    else
        linespec = [0 0 0];
    end

    plot(burstIndex(useCC & rgn),burstISI(useCC & rgn),'o','MarkerEdgeColor',linespec)
end
ylim([0 8]);
ylabel('Burst ISI (ms)')
xscale log; xlim([0 20]); grid on
xlabel('Burst Index')
set(gca,'FontSize',12,'FontName','Arial')

% fhandle = figure;
% set(gcf,'units','normalized','position',[0.4 0.35 0.20 0.39])
% fhandle = fixRatio(fhandle);
% plot3(burstIndex(useCC),burstISI(useCC),burstLens(useCC),'ko')
% ylim([0 8]);
% xlabel('Burst Index (Burst spikes / Reg spikes)')
% ylabel('Mean Burst ISI (ms)')
% zlabel('Mean Burst length')
% set(gca,'FontSize',12,'FontName','Arial')

end