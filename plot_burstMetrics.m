function [fhandle] = plot_burstMetrics(root,sess,useCC,rgdat,rgncols)

arguments
    root
    sess
    useCC
    rgdat
    rgncols = [0.9961 0.7305 0.4336; 0.3672 0.2969 0.3711; 0 0 0]; % Sub gold vs CA1 dull purple
end

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
        linespec = rgncols(1,:);
    elseif nrgns(i) == 2    % Sub
        linespec = rgncols(2,:);
    else
        linespec = rgncols(3,:);
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