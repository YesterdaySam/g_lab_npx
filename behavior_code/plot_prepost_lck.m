function [fhandle] = plot_prepost_lck(sess1,sess2,dbnsz,plotflag)
%% Overlay linearized velocity (binned by space) pre and post
% Inputs
%   sess        = struct from importBhvr.m
%   dbnsz        = double in meters (m)
% Outputs
%   fhandle     = handle to figure 1

arguments
    sess1
    sess2
    dbnsz = 0.03    % velocity bin size in m
    plotflag = 1
end

[binedges1,~,pslck1] = plot_lickpos(sess1,dbnsz,0);
[binedges2,~,pslck2] = plot_lickpos(sess2,dbnsz,0);
normlck1 = pslck1 ./ sum(pslck1,2);
normlck2 = pslck2 ./ sum(pslck2,2);
r1pos = 100*round(mean(sess1.pos(sess1.rwdind(1:30))),1);
r2pos = 100*round(mean(sess2.pos(sess2.rwdind(1:30))),1);

if plotflag
    tmp1 = figure;
    set(gcf,'units','normalized','position',[0.4 0.35 0.25 0.45])
    imagesc(pslck1,[prctile(pslck1,1,'all'), prctile(pslck1,99,'all')]);
    colormap("turbo")
    cbar = colorbar; clim([0 inf]);
    xlabel('Position'); % xlim([0 200])
    xticks(1:20:length(binedges1)); xticklabels(binedges1(1:20:length(binedges1))*100);
    ylabel('Trial #'); ylabel(cbar,'Licks/s','FontSize',12,'Rotation',90)
    set(gca,'FontSize',12,'FontName','Arial')

    tmp2 = figure;
    set(gcf,'units','normalized','position',[0.4 0.35 0.25 0.45])
    imagesc(pslck2,[prctile(pslck2,1,'all'), prctile(pslck2,99,'all')]);
    colormap("turbo")
    cbar = colorbar; clim([0 inf]);
    xlabel('Position'); % xlim([0 200])
    ylabel('Trial #'); ylabel(cbar,'Licks/s','FontSize',12,'Rotation',90)
    xticks(1:20:length(binedges2)); xticklabels(binedges2(1:20:length(binedges2))*100);
    set(gca,'FontSize',12,'FontName','Arial')

    sem1 = std(normlck1,'omitnan')/sqrt(sess1.nlaps);
    ciup1 = rmmissing(mean(normlck1,1,'omitnan') + sem1*1.96);
    cidn1 = rmmissing(mean(normlck1,1,'omitnan') - sem1*1.96);
    sem2 = std(normlck2,'omitnan')/sqrt(sess2.nlaps);
    ciup2 = rmmissing(mean(normlck2,1,'omitnan') + sem2*1.96);
    cidn2 = rmmissing(mean(normlck2,1,'omitnan') - sem2*1.96);

    fhandle = figure; hold on
    ylim([0 0.2])
    set(gcf,'units','normalized','position',[0.4 0.35 0.22 0.35])
    patch(100*[binedges1(1:length(cidn1)),fliplr(binedges1(1:length(cidn1)))],[cidn1,fliplr(ciup1)],...
        'k','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(binedges1(1:end-1)*100,mean(normlck1,1,'omitnan'),'Color','k','LineWidth',2)
    plot([r1pos, r1pos], ylim,'k--','HandleVisibility','off');

    patch(100*[binedges2(1:length(cidn2)),fliplr(binedges2(1:length(cidn2)))],[cidn2,fliplr(ciup2)],...
        'r','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(binedges2(1:end-1)*100,mean(normlck2,1,'omitnan'),'Color','r','LineWidth',2)
    plot([r2pos, r2pos], ylim,'r--','HandleVisibility','off');
    
    xlabel('Position'); xlim([0 100*max(binedges1)])
    ylabel('Lick probability')
    legend('Familiar RZ','Novel RZ')
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end
end