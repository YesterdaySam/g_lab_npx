% Firing raster for neural only recordings
ccs = root.info.cluster_id(root.goodind & root.info.fr < 20 & root.info.uType);

% Seconds
winstt = 5;
winend = 35;

for i = 1:length(ccs)
    spkstruc(i).spks = root.tsb(root.cl == ccs(i));
    spkstruc(i).lapspks = spkstruc(i).spks(spkstruc(i).spks > winstt & spkstruc(i).spks < winend);
end

%%

rastf = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
% plot(sess.ts(exlapstt:exlapend),sess.pos(exlapstt:exlapend)./pmax,'b')
% plot(sess.ts(exlapstt:exlapend),1.15+sess.velshft(exlapstt:exlapend)./vmax,'r')
cmap = jet(length(ccs));

for i = 1:length(ccs)
    plot(spkstruc(i).lapspks - winstt, ones(length(spkstruc(i).lapspks),1) - 1 + 1*i, '|','Color',cmap(i,:))
end

xlabel('Time (s)')
ylabel("Neuron #")
ylim([0 inf])
% legend({'Position','Velocity','Spikes'})
set(gca,'FontSize',12,'FontName','Arial')

fsave(rastf,[root.name '_goodRaster'],1,1)