function [binedges,bnvel,fhandle,fhandle2] = plot_trialvel(sess,dbnsz,plotflag)
%% Create linearized velocity (binned by space)
% Inputs
%   sess        = struct from importBhvr.m
%   bnsz        = double in meters (m)
% Outputs
%   fhandle     = handle to figure 1
%   fhandle2    = handle to figure 2
%   trackedges  = velocity bin edges
%   bnvel       = binned velocity

arguments
    sess
    dbnsz = 0.01    % velocity bin size in m
    plotflag = 1
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

try 
    binedges = 0:dbnsz:sess.maxPos;
catch
    binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
end
nbins           = size(binedges,2);

for i = 1:sess.nlaps
    bnocc = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges);
    for j = 1:nbins-1
        tmpind = sess.lapstt(i) + find(sess.pos(sess.lapstt(i):sess.lapend(i)) > binedges(j) & sess.pos(sess.lapstt(i):sess.lapend(i)) < binedges(j+1));
        bnvel(i,j) = mean(sess.velshft(tmpind(1:end-1)));    % Average everything except the last time index which bleeds into next lap (or over recording length)
    end
    % bnvel(i,:) = bnvel(i,:) ./ (bnocc / sess.samprate);
end

if plotflag
    fhandle = figure;
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    imagesc(bnvel,[prctile(bnvel,1,'all'), prctile(bnvel,99,'all')]);
    % imagesc(bnvel, [0 max(bnvel,[],'all')])
    colormap("winter")
    cbar = colorbar; clim([0 inf]);
    xlabel('Position'); % xlim([0 200])
    % xticks(1:45:length(binedges)); xticklabels(binedges(1:45:length(binedges))*100);
    ylabel('Trial #'); ylabel(cbar,'cm/s','FontSize',12,'Rotation',90)
    set(gca,'FontSize',12,'FontName','Arial','YDir','reverse')

    sem = std(bnvel,'omitnan')/sqrt(sess.nlaps);
    ciup = rmmissing(mean(bnvel,1,'omitnan') + sem*1.96);
    cidn = rmmissing(mean(bnvel,1,'omitnan') - sem*1.96);

    fhandle2 = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
    plot(binedges(1:end-1)*100,mean(bnvel,1,'omitnan'),'k','LineWidth',2)
    patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    % plot(bnvel','Color',[.5 .5 .5])
    xlabel('Position'); % xlim([0 200])
    ylabel('Average Velocity'); ylim([0 prctile(sess.velshft,99)])
    set(gca,'FontSize',12,'FontName','Arial','YDir','reverse')
end

end