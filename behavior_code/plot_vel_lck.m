function [fhandle] = plot_vel_lck(sess,bnvel,bnlck,edges1,edges2)

nLaps = size(bnvel,1);   %Recalculate basedo on laps used in bnvel, don't use sess.nlaps

semvel = std(bnvel,'omitnan')/sqrt(nLaps);
ciupvel = rmmissing(mean(bnvel,1,'omitnan') + semvel*1.96);
cidnvel = rmmissing(mean(bnvel,1,'omitnan') - semvel*1.96);

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
plot(edges1(1:end-1)*100,mean(bnvel,1,'omitnan'),'k','LineWidth',2)
patch(100*[edges1(1:length(cidnvel)),fliplr(edges1(1:length(cidnvel)))],[cidnvel,fliplr(ciupvel)],'k','FaceAlpha',0.5,'EdgeColor','none')
xlabel('Position'); xlim([0 200])
ylabel('Average Velocity'); ylim([0 prctile(sess.velshft,99)*100])
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

semlck = std(bnlck,'omitnan')/sqrt(nLaps);
ciuplck = rmmissing(mean(bnlck,1,'omitnan') + semlck*1.96);
cidnlck = rmmissing(mean(bnlck,1,'omitnan') - semlck*1.96);

yyaxis right
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
plot(edges2(1:end-1)*100,mean(bnlck,1,'omitnan'),'b','LineWidth',2)
patch(100*[edges2(1:length(cidnlck)),fliplr(edges2(1:length(cidnlck)))],[cidnlck,fliplr(ciuplck)],'b','FaceAlpha',0.5,'EdgeColor','none')
ylabel('Average Licks/s'); ylim([0 inf])
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

tmp = gca;
tmp.YAxis(2).Color = 'b';
tmp.YAxis(1).Color = 'k';

end