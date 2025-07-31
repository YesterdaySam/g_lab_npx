function [fhandle] = plot_prepost_vel(sess1,sess2,dbnsz,plotflag)
%% Overlay linearized velocity (binned by space) pre and post
% Inputs
%   sess        = struct from importBhvr.m
%   dbnsz        = double in meters (m)
% Outputs
%   fhandle     = handle to figure 1

arguments
    sess1
    sess2
    dbnsz = 0.01    % velocity bin size in m
    plotflag = 1
end

[binedges1,bnvel1] = plot_trialvel(sess1,dbnsz,0);
[binedges2,bnvel2] = plot_trialvel(sess2,dbnsz,0);
r1pos = 100*round(mean(sess1.pos(sess1.rwdind(1:30))),1);
r2pos = 100*round(mean(sess2.pos(sess2.rwdind(1:30))),1);

if plotflag
    sem1 = std(bnvel1,'omitnan')/sqrt(sess1.nlaps);
    ciup1 = rmmissing(mean(bnvel1,1,'omitnan') + sem1*1.96);
    cidn1 = rmmissing(mean(bnvel1,1,'omitnan') - sem1*1.96);
    sem2 = std(bnvel2,'omitnan')/sqrt(sess2.nlaps);
    ciup2 = rmmissing(mean(bnvel2,1,'omitnan') + sem2*1.96);
    cidn2 = rmmissing(mean(bnvel2,1,'omitnan') - sem2*1.96);

    fhandle = figure; hold on
    ylim([0 prctile(sess1.velshft,99)])
    cmapwinter = winter(2);
    set(gcf,'units','normalized','position',[0.4 0.35 0.22 0.35])
    patch(100*[binedges1(1:length(cidn1)),fliplr(binedges1(1:length(cidn1)))],[cidn1,fliplr(ciup1)],...
        'b','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(binedges1(1:end-1)*100,mean(bnvel1,1,'omitnan'),'Color',cmapwinter(1,:),'LineWidth',2)
    plot([r1pos, r1pos], ylim,'b--','HandleVisibility','off');

    patch(100*[binedges2(1:length(cidn2)),fliplr(binedges2(1:length(cidn2)))],[cidn2,fliplr(ciup2)],...
        'g','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(binedges2(1:end-1)*100,mean(bnvel2,1,'omitnan'),'Color',cmapwinter(2,:),'LineWidth',2)
    plot([r2pos, r2pos], ylim,'g--','HandleVisibility','off');
    
    xlabel('Position'); xlim([0 100*max(binedges2)])
    ylabel('Velocity (cm/s)')
    legend('Familiar RZ','Novel RZ')
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end
end