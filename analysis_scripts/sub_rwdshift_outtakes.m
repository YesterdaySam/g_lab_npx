%% Averaged lick and velocity profiles

for i = 1:length(mID)
    uPreVel(i,:) = mean(bvDat(i).preVMap,'omitnan');
    uPstVel(i,:) = mean(bvDat(i).pstVMap,'omitnan');
    uPreVF20(i,:) = bvDat(i).uPreVF20;
    uPreVL20(i,:) = bvDat(i).uPreVL20;
    uPstVF20(i,:) = bvDat(i).uPstVF20;
    uPstVL20(i,:) = bvDat(i).uPstVL20;

    uPreLck(i,:) = mean(bvDat(i).preLMap,'omitnan');
    uPstLck(i,:) = mean(bvDat(i).pstLMap,'omitnan');
    uPreLF20(i,:) = bvDat(i).uPreLF20;
    uPreLL20(i,:) = bvDat(i).uPreLL20;
    uPstLF20(i,:) = bvDat(i).uPstLF20;
    uPstLL20(i,:) = bvDat(i).uPstLL20;
end

uPreVelF = plot_bhvrTraceXmouse(circshift(uPreVel,30,2),uPreVF20,uPstVL20,[0.45 0.45 0.45; 0 0 0]);
xticks([0 40 80 120 160]); xticklabels([-30 10 50 90 130]);
ylabel('Velocity (cm/s)'); ylim([0 50]); xlim([0 186])
patch([xcoordsPre,fliplr(xcoordsPre)],[0 0 50 50],[0 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
patch([xcoordsPst,fliplr(xcoordsPst)],[0 0 50 50],[1 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
plot([40 40],[0 50],'--','Color',vColors2(1,:)); plot([130 130],[0 50],'--','Color',vColors2(2,:))

uPstVelF = plot_bhvrTraceXmouse(circshift(uPstVel,30,2),uPstVF20,uPstVL20,[1 .65 .65; 1 0 0]);
xticks([0 40 80 120 160]); xticklabels([-30 10 50 90 130]);
ylabel('Velocity (cm/s)'); ylim([0 50]); xlim([0 186])
patch([xcoordsPre,fliplr(xcoordsPre)],[0 0 50 50],[0 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
patch([xcoordsPst,fliplr(xcoordsPst)],[0 0 50 50],[1 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
plot([40 40],[0 50],'--','Color',vColors2(1,:)); plot([130 130],[0 50],'--','Color',vColors2(2,:))

uPreLckF = plot_bhvrTraceXmouse(circshift(uPreLck,10,2),uPreLF20,uPreLL20,[0.45 0.45 0.45; 0 0 0]);
xticks([0 40 80 120 160]); xticklabels([-30 10 50 90 130]);
ylabel('Licks/spatial bin'); ylim([0 12]); xlim([0 186])
patch([xcoordsPre,fliplr(xcoordsPre)],[0 0 12 12],[0 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
patch([xcoordsPst,fliplr(xcoordsPst)],[0 0 12 12],[1 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
plot([40 40],[0 12],'--','Color',vColors2(1,:)); plot([130 130],[0 12],'--','Color',vColors2(2,:))

uPstLckF = plot_bhvrTraceXmouse(circshift(uPstLck,10,2),uPstLF20,uPstLL20,[1 .65 .65; 1 0 0]);
xticks([0 40 80 120 160]); xticklabels([-30 10 50 90 130]);
ylabel('Licks/spatial bin'); ylim([0 12]); xlim([0 186])
patch([xcoordsPre,fliplr(xcoordsPre)],[0 0 12 12],[0 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
patch([xcoordsPst,fliplr(xcoordsPst)],[0 0 12 12],[1 0 0],'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
plot([40 40],[0 12],'--','Color',vColors2(1,:)); plot([130 130],[0 12],'--','Color',vColors2(2,:))

if saveFlag
    fsave(uPreVelF,[sbase 'bhv_meanVelPre'])
    fsave(uPstVelF,[sbase 'bhv_meanVelPst'])
    fsave(uPreLckF,[sbase 'bhv_meanLckPre'])
    fsave(uPstLckF,[sbase 'bhv_meanLckPst'])
end

%% Group analyses

%% PVC Sparseness
% Sparseness - P(alpha) From Battaglia et al., 2004
% Normalize all PF firing rates, sum by bin, and sqaure
% Divide by sum of squared PF firing rates by bin, and normalize by N cells

numerPre = sum(posNormPre).^2;
denomPre = sum(posNormPre.^2);

pvSparsePre = numerPre ./ denomPre ./ size(posNormPre,1);



%%
% function [fhandle] = plotDistroHisto(distro1,distro2,binpos,rzPos)
% vColors2 = [0.5 0.5 1; 0.75 0.75 1];
% 
% fhandle = figure; hold on
% plot(binpos*100,distro1./sum(distro1),'Color',vColors2(1,:),'LineWidth',2);
% plot(binpos*100,distro2./sum(distro2),'Color',vColors2(2,:),'LineWidth',2);
% if length(rzPos) == 2
%     plot([rzPos(1) rzPos(1)]*100,[0 0.15],'--','Color',vColors2(1,:))
%     plot([rzPos(2) rzPos(2)]*100,[0 0.15],'--','Color',vColors2(2,:))
% else
%     plot([rzPos rzPos]*100,[0 0.15],'k--')
% end
% ylim([0 max([distro1./sum(distro1); distro2./sum(distro2)+0.02],[],'all')])
% xlim([binpos(1)*100-1 binpos(end)*100+1])
% ylabel('Probability of field peak')
% legend({'Familiar','Novel'},'Location','northeast')
% set(gca,'FontSize',12,'FontName','Arial')
% end

% function [fhandle] = plot_barXmouse(dat)
% vColors2 = [0.5 0.5 1; 0.75 0.75 1];
% 
% fhandle = figure; hold on
% set(gcf,'units','normalized','position',[0.4 0.35 0.12 0.39])
% b = bar(mean(dat),'FaceColor','flat','HandleVisibility','off');
% b.CData = vColors2;
% errorbar(1:2,mean(dat),std(dat)/sqrt(size(dat,1)),'k.','HandleVisibility','off')
% plot(dat','k-o')
% xlim([0.5 2.5])
% xticks(1:2); xticklabels({'Familiar', 'Novel'})
% legend('Mouse')
% set(gca,'FontSize',12,'FontName','Arial')
% end

function [fhandle] = plot_bhvrTraceXmouse(dat,datF20,datL20,vcolors)

bnpos = linspace(0,185,size(dat,2));

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.22 0.4])
plot(bnpos,dat','Color',vcolors(1,:))
plot(bnpos,mean(dat),'Color',vcolors(2,:),'LineWidth',2)
% plot(mean(datF20),'Color',vcolors(1,:),'LineWidth',2);
% plot(mean(datL20),'Color',vcolors(2,:),'LineWidth',2);
xlabel('Position (cm)')
set(gca,'FontSize',16,'FontName','Arial')

end

