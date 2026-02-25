function [binedges,bnvel,fhandle] = get_velXpos(sess,dbnsz,plotflag)
%% Create linearized velocity (binned by time)
% !! Not functional!
% Inputs
%   sess        = struct from importBhvr.m
%   dbnsz       = double in meters (m)
%   plotFlag    = binary, whether to plot
% Outputs
%   binedges    = spatial bin edges
%   bnvel       = avg velocity by spatial bin
%   fhandle     = figure
%
% Created 7/7/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    dbnsz = 0.01    % velocity bin size in m
    plotflag = 1
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
nbins           = size(binedges,2);
clear bnvel

% for i = 1:sess.nlaps
%     tmpVel(i,:) = sess.velshft(binedges(1)+sess.rwdind(i):binedges(end)+sess.rwdind(i));
% end

for i = 1:length(binedges)-1
    posIncl = sess.pos >= binedges(i) & sess.pos < binedges(i+1);
    bnvel(i) = mean(sess.velshft(posIncl));
end
bnocc = histcounts(sess.pos(sess.lapstt(1):sess.lapend(end)),binedges);

if plotflag
    fhandle = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
    plot(binedges(1:end-1)*100,bnvel,'k','LineWidth',2)
    % patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    % plot(bnvel','Color',[.5 .5 .5])
    xlabel('Position'); % xlim([0 200])
    ylabel('Average Velocity'); ylim([0 prctile(sess.velshft,99)])
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end