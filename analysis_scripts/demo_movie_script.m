%% Make demo video of behavior, video, and some spiking cells

load('D:\Data\Kelton\analyses\KW022\KW022_12172024_rec_D2_RMed1\KW022_12172024_rec_D2_RMed1_root.mat') % Root file
load('D:\Data\Kelton\analyses\KW022\KW022_12172024_rec_D2_RMed1\KW022_12172024_session.mat') % Sess file
load('D:\Data\Kelton\analyses\KW022\KW022_12172024_rec_D2_RMed1\KW022_12172024_rec_D2_RMed1_goodPFs.mat') % GoodPfs file
vReader = VideoReader("D:\Data\Kelton\bhvr_data\KW022\KW022_12172024_rec_D2_RMed1_bhvr.avi");

%% Test frame of behavior

% Test frame of behavior
mLen = floor(vReader.Duration);
maxPos = max(sess.pos);
maxVel = max(sess.velshft);
vSmth  = smoothdata(sess.velshft,'Gaussian',1000);
nInds = sess.samprate*mLen;

% goodshank = 2;
% ca1PFs = goodPFs(:,3) == goodshank;
% ca1PFs = goodPFs(ca1PFs,:);
units = [1 13 60 110 111 324 477 519 520];
jetmap = jet(length(units));

figure; hold on; axis off
set(gcf,'color','w','units','normalized','position',[0.2 0.1 0.5 0.8])
% set(gcf,'color','w','units','normalized','position',[0.2 0.1 0.35 0.3])
plot(sess.ts(1:nInds), 10*(vSmth(1:nInds)./maxVel),'r')
plot(sess.ts(1:nInds), 10*(sess.pos(1:nInds)./maxPos)+10,'b')
tmplcks = sess.lckind(sess.lckind < nInds);
plot(sess.ts(tmplcks),10 + zeros(length(tmplcks),1),'k|')
for i = 1:length(units) 
    tmpspks = root.cl == units(i);
    tmpspks = root.tsb(tmpspks);
    dispspks = tmpspks(tmpspks < nInds);
    plot(sess.ts(dispspks),(20+i*3) + ones(1,length(dispspks))','|','Color',jetmap(i,:))
end
ylim([-0.5 190])
% xlabel('Time')
% ylabel('Normalized value')
% legend({'Velocity','Position'})
set(gca,'FontSize',12,'FontName','Arial')

%% Plot a few hand-selected neurons relative to behavior

figure; hold on; axis off
set(gcf,'color','w','units','normalized','position',[0.2 0.1 0.5 0.3])
xlim([0 sess.ts(nInds)])
ylim([-.2 4])
plot(sess.ts(1:10000), vSmth(1:10000)./maxVel,'r')
plot(sess.ts(1:10000), sess.pos(1:10000)./maxPos+1.15,'b')
for i = 1:length(units)
    tmpspks = root.cl == units(i);
    tmpspks = root.tsb(tmpspks);
    dispspks = tmpspks(tmpspks < 10000);
    plot(sess.ts(dispspks),(1.2+i/10) + ones(1,length(dispspks))','|','Color',jetmap(i,:))
end
% legend({'Velocity','Position','Spike Raster'})
set(gca,'FontSize',12,'FontName','Arial')

set(gca,"NextPlot","replacechildren")

v = VideoWriter("test.avi");
open(v)

vStart = round(sess.ts(sess.lapstt(6))*30);
iStart = sess.lapstt(6);
vEnd = round(sess.ts(sess.lapend(9))*30);
iEnd = sess.lapend(9);
for k = vStart:vEnd
    subFrame = round(k / (30/sess.samprate));
    plot(sess.ts(iStart:subFrame), vSmth(iStart:subFrame)./maxVel,'r')
    xlim([sess.ts(iStart) sess.ts(iEnd)])
    ylim([-.2 4])

    hold on
    plot(sess.ts(iStart:subFrame), sess.pos(iStart:subFrame)./maxPos+1.15,'b')
    tmplcks = sess.lckind(sess.lckind > iStart & sess.lckind < subFrame);
    plot(sess.ts(tmplcks), 0.8 + zeros(length(tmplcks),1),'k|')

    for i = 1:length(units)
        tmpspks = root.cl == units(i);
        tmpspks = root.tsb(tmpspks);
        dispspks = tmpspks(tmpspks > iStart & tmpspks < subFrame);
        plot(sess.ts(dispspks),(1.2+i/5) + ones(1,length(dispspks))','|','Color',jetmap(i,:))
    end
    % legend({'Velocity','Position','Spike Raster'})
    set(gca,'FontSize',12,'FontName','Arial')
    axis off
    frame = getframe(gcf);
    writeVideo(v,frame)
    hold off

end
close(v)

%% Plot many neurons relative to Position only
% [~, tmpSortInds] = plot_unitsXpos(root,sess,root.good);
hsvmap = hsv(length(root.good));

figure; hold on; axis off
set(gcf,'color','w','units','normalized','position',[0.2 0.1 0.5 0.8])

% xlim([sess.ts(iStart) sess.ts(iEnd)])
% ylim([0 length(root.good) + 5])
% plot(sess.ts(iStart:iEnd), 5*(sess.pos(iStart:iEnd)./maxPos),'b')
% for i = 1:length(root.good)
%     tmpspks = root.cl == root.good(find(tmpSortInds == i));
%     tmpspks = root.tsb(tmpspks);
%     dispspks = tmpspks(tmpspks > iStart & tmpspks < iEnd);
%     plot(sess.ts(dispspks), 5 + i + ones(1,length(dispspks))','|','Color',hsvmap(i,:))
% end
% % legend({'Velocity','Position','Spike Raster'})
% set(gca,'FontSize',12,'FontName','Arial')

set(gca,"NextPlot","replacechildren")

v = VideoWriter("demo_GoodUnits_lap6-9.avi");
open(v)

lStt = 6;
lEnd = 9;
vStart = round(sess.ts(sess.lapstt(lStt))*30);
iStart = sess.lapstt(lStt);
vEnd = round(sess.ts(sess.lapend(lEnd))*30);
iEnd = sess.lapend(lEnd);
for k = vStart:vEnd
    subFrame = round(k / (30/sess.samprate));
    xlim([sess.ts(iStart) sess.ts(iEnd)])
    ylim([0 length(root.good) + 5])

    hold on
    plot(sess.ts(iStart:subFrame), 5*(sess.pos(iStart:subFrame)./maxPos),'b')

    for i = 1:length(root.good)
        tmpspks = root.cl == root.good(find(tmpSortInds == i));
        tmpspks = root.tsb(tmpspks);
        dispspks = tmpspks(tmpspks > iStart & tmpspks < subFrame);
        plot(sess.ts(dispspks), 5 + i + ones(1,length(dispspks))','|','Color',hsvmap(i,:))
    end
    set(gca,'FontSize',12,'FontName','Arial')
    axis off
    frame = getframe(gcf);
    writeVideo(v,frame)
    hold off

end
close(v)

%% Plot many neurons relative to behavior

hsvmap = hsv(length(root.good));

figure; hold on; axis off
set(gcf,'color','w','units','normalized','position',[0.2 0.1 0.5 0.8])
set(gca,"NextPlot","replacechildren")
ylim([-1 190])

v = VideoWriter("KW022_demo_laps6-12.avi");
open(v)

lStt = 6;
lEnd = 9;
vStart = round(sess.ts(sess.lapstt(lStt))*30);
iStart = sess.lapstt(lStt);
vEnd = round(sess.ts(sess.lapend(lEnd))*30);
iEnd = sess.lapend(lEnd);
for k = vStart:vEnd
    subFrame = round(k / (30/sess.samprate));
    plot(sess.ts(iStart:subFrame), 10*vSmth(iStart:subFrame)./maxVel,'r')
    xlim([sess.ts(iStart) sess.ts(iEnd)])
    ylim([-1 190])

    hold on
    plot(sess.ts(iStart:subFrame), 10*(sess.pos(iStart:subFrame)./maxPos)+10,'b')
    tmplcks = sess.lckind(sess.lckind > iStart & sess.lckind < subFrame);
    plot(sess.ts(tmplcks), 10 + zeros(length(tmplcks),1),'k|')

    for i = 1:length(units)
        tmpspks = root.tsb(root.cl == units(i));
        tmpspks = root.tsb(tmpspks);
        dispspks = tmpspks(tmpspks > iStart & tmpspks < subFrame);
        plot(sess.ts(dispspks),(20+i*3) + ones(1,length(dispspks))','|','Color',jetmap(i,:))
    end
    % legend({'Velocity','Position','Spike Raster'})
    set(gca,'FontSize',12,'FontName','Arial')
    axis off
    frame = getframe(gcf);
    writeVideo(v,frame)
    hold off
end

lStt = 10;
lEnd = 12;
vStart = round(sess.ts(sess.lapstt(lStt))*30);
iStart = sess.lapstt(lStt);
vEnd = round(sess.ts(sess.lapend(lEnd))*30);
iEnd = sess.lapend(lEnd);

for k = vStart:vEnd
    subFrame = round(k / (30/sess.samprate));
    plot(sess.ts(iStart:subFrame), 10*vSmth(iStart:subFrame)./maxVel,'r')
    xlim([sess.ts(iStart) sess.ts(iEnd)])
    ylim([-1 190])

    hold on
    plot(sess.ts(iStart:subFrame), 10*(sess.pos(iStart:subFrame)./maxPos)+10,'b')
    tmplcks = sess.lckind(sess.lckind > iStart & sess.lckind < subFrame);
    plot(sess.ts(tmplcks), 10 + zeros(length(tmplcks),1),'k|')

    for i = 1:length(root.good)
        tmpspks = root.tsb(root.cl == root.good(i));
        dispspks = tmpspks(tmpspks > iStart & tmpspks < subFrame);
        plot(sess.ts(dispspks), 20 + i + ones(1,length(dispspks))','|','Color',hsvmap(i,:))
    end
    set(gca,'FontSize',12,'FontName','Arial')
    axis off
    frame = getframe(gcf);
    writeVideo(v,frame)
    hold off
end
close(v)

%% Old scrap code -- just combine videos with MS clipchamp as post production
% % Test frame of video
% s = struct("cdata",zeros(vidObj.Height,vidObj.Width,3,"uint8"),colormap=[]);
% 
% vReader.CurrentTime = 30;
% vFrame = readFrame(vReader);
% figure; 	axis('tight')
% set(gcf,'units','normalized','position',[0.2 0.1 0.7 0.7])
% imshow(vFrame)
% 
% % Merge together in subplots
% figlist = get(groot, 'Children');
% newfig = figure;
% set(gcf,'units','normalized','position',[0.2 0.1 0.4 0.7])
% 
% tcl = tiledlayout(newfig, 'flow');
% 
% for j = 1:numel(figlist)
%     figure(figlist(j))
%     ax = gca;
%     ax.Parent = tcl;
%     ax.Layout.Tile = j;
% end
% title(tcl,['Unit ' num2str(cc)])
% 


%% 
Z = peaks;
surf(Z);
axis tight manual
set(gca,"NextPlot","replacechildren")
v = VideoWriter("test.avi");
open(v)

for k = 1:20
   surf(sin(2*pi*k/20)*Z,Z)
   frame = getframe(gcf);
   writeVideo(v,frame)
end

close(v)