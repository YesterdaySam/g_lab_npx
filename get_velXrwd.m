function [binedges,bnvel,fhandle] = get_velXrwd(sess,tbnsz,wlen,plotflag)
%% Create linearized velocity (binned by time)
% Inputs
%   sess        = struct from importBhvr.m
%   tbnsz       = double in seconds (s)
%   wlen        = double in seconds (s)
%   plotFlag    = binary, whether to plot
% Outputs
%   binedges    = time bin edges
%   bnvel       = avg velocity by time bin
%   fhandle     = figure
%
% Created 7/7/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    tbnsz = 0.25    % time bin size in s
    wlen  = 5
    plotflag = 1
end

tbnInd = tbnsz / (1 / sess.samprate);
wlenInd = wlen / (1 / sess.samprate);
binedges = -wlenInd:tbnInd:wlenInd;

for i = 1:sess.nlaps
    try
        allvel(i,:) = sess.velshft(binedges(1)+sess.rwdind(i):binedges(end)+sess.rwdind(i));
    catch
        % disp(['Skipping lap ' num2str(i)]);
    end
end
% bnocc = histcounts(sess.pos(sess.lapstt(1):sess.lapend(end)),binedges);
[~,~,loc] = histcounts(-wlenInd:wlenInd,binedges);
bnvel = accumarray(loc(:),mean(allvel)) ./ accumarray(loc(:),1);

if plotflag
    sem = rmmissing(std(mean(allvel),'omitnan')/sqrt(sess.nlaps));
    ciup = bnvel + sem*1.96;
    cidn = bnvel - sem*1.96;

    fhandle = figure; hold on
    xcoords = (binedges(1:end-1)+0.5*tbnInd)/sess.samprate;
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
    plot(xcoords,bnvel,'k','LineWidth',2)
    patch([xcoords,fliplr(xcoords)],[cidn',fliplr(ciup')],'k','FaceAlpha',0.5,'EdgeColor','none')
    xlabel('Time to reward (s)');
    ylabel('Average Velocity'); ylim([0 prctile(sess.velshft,99)])
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end