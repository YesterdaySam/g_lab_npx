% Plot example firing raster over laps

% ccs(1) = 160;
% ccs(2) = 154;
% ccs(3) = 143;
ccs = root.good;

exlapstt = sess.lapstt(22);
exlapend = sess.lapend(23);

for i = 1:length(ccs)
    spkstruc(i).spks = sess.ts(root.tsb(root.cl == ccs(i)));
    spkstruc(i).lapspks = spkstruc(i).spks(spkstruc(i).spks > sess.ts(exlapstt) & spkstruc(i).spks < sess.ts(exlapend));
end

vmax = max(sess.velshft(exlapstt:exlapend));
pmax = max(sess.pos(exlapstt:exlapend));


%%

figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.5 0.3])
% plot(sess.ts(exlapstt:exlapend),sess.pos(exlapstt:exlapend)./pmax,'b')
% plot(sess.ts(exlapstt:exlapend),1.15+sess.velshft(exlapstt:exlapend)./vmax,'r')
cmap = jet(length(root.good));

for i = 1:length(ccs)
    plot(spkstruc(i).lapspks, ones(length(spkstruc(i).lapspks),1) + 1 + 1*i, '|','Color',cmap(i,:))
end

xlabel('Time (s)')
ylabel("Neuron #")
ylim([0 inf])
% legend({'Position','Velocity','Spikes'})
set(gca,'FontSize',12,'FontName','Arial')