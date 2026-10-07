function [preHeatmapF,pstHeatmapF,averageSumF] = plot_prepost_lck(sess1,sess2,dbnsz,plotflag)
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
pslck1mask = pslck1 > 0;
pslck2mask = pslck2 > 0;
normlck1 = sum(pslck1mask) ./ size(pslck1,1); % Normalizing by P(lck) over laps and spatial bins
normlck2 = sum(pslck2mask) ./ size(pslck2,1); 
% normlck1 = pslck1 ./ sum(pslck1,2); % Normalizing by space
% normlck2 = pslck2 ./ sum(pslck2,2);

try
    r1pos = 100*round(mean(sess1.pos(sess1.rwdind(1:30))),1);
catch
    disp('Fewer than 30 trials before shift, using all')
    r1pos = 100*round(mean(sess1.pos(sess1.rwdind(1:end))),1);
end
try
    r2pos = 100*round(mean(sess2.pos(sess2.rwdind(1:30))),1);
catch
    disp('Fewer than 30 trials after shift, using all')
    r2pos = 100*round(mean(sess2.pos(sess2.rwdind(1:end))),1);
end

if plotflag
    preHeatmapF = plot_lckheatmap(sess1);
    colorbar off

    pstHeatmapF = plot_lckheatmap(sess2);
    colorbar off

    sem1 = std(normlck1,'omitnan')/sqrt(sess1.nlaps);
    ciup1 = rmmissing(mean(normlck1,1,'omitnan') + sem1*1.96);
    cidn1 = rmmissing(mean(normlck1,1,'omitnan') - sem1*1.96);
    sem2 = std(normlck2,'omitnan')/sqrt(sess2.nlaps);
    ciup2 = rmmissing(mean(normlck2,1,'omitnan') + sem2*1.96);
    cidn2 = rmmissing(mean(normlck2,1,'omitnan') - sem2*1.96);

    averageSumF = figure; hold on
    ylim([0 1])
    set(gcf,'units','normalized','position',[0.4 0.35 0.22 0.35])
    fixRatio(averageSumF);
    patch(100*[binedges1(1:length(cidn1)),fliplr(binedges1(1:length(cidn1)))],[cidn1,fliplr(ciup1)],...
        'k','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(binedges1(1:end-1)*100,mean(normlck1,1,'omitnan'),'Color','k','LineWidth',2)
    plot([r1pos, r1pos], ylim,'k--','HandleVisibility','off');

    patch(100*[binedges2(1:length(cidn2)),fliplr(binedges2(1:length(cidn2)))],[cidn2,fliplr(ciup2)],...
        'r','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(binedges2(1:end-1)*100,mean(normlck2,1,'omitnan'),'Color','r','LineWidth',2)
    plot([r2pos, r2pos], ylim,'r--','HandleVisibility','off');
    
    xlabel('Position'); xlim([0 100*max(binedges1)])
    ylabel('P(Lick x spatial bin) across laps')
    legend('Familiar RZ','Novel RZ')
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end
end