function [fhandle,fhandle2, trackedges, bnvel] = plot_trialvel(sess, bnsz)
%% Create linearized velocity (binned by space)
% Inputs
% session = struct from importBhvr.m
% bnsz    = double in cm e.g. 1
% Outputs
% fhandle = handle to figure

arguments
    sess
    bnsz = 1    % cm
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

% vel             = session.vel - min(session.vel);
vel             = sess.velshft;
bnsz            = bnsz/100;    % translate to m
tracklen        = (max(sess.pos) - min(sess.pos)); % m
trackedges      = 0:bnsz:tracklen;
nbins           = size(trackedges,2);
clear bnvel

for i = 1:sess.nlaps
    bnocc = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),trackedges);
    for j = 1:nbins-1
        tmpind = sess.lapstt(i) + find(sess.pos(sess.lapstt(i):sess.lapend(i)) > trackedges(j) & sess.pos(sess.lapstt(i):sess.lapend(i)) < trackedges(j+1));
        bnvel(i,j) = mean(vel(tmpind))*100;
    end
    % bnvel(i,:) = bnvel(i,:) ./ (bnocc / sess.samprate);
end

fhandle = figure;
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
imagesc(bnvel,[prctile(bnvel,1,'all'), prctile(bnvel,99,'all')]); 
% imagesc(bnvel, [0 max(bnvel,[],'all')])
colormap("sky")
cbar = colorbar; clim([0 inf]);
xlabel('Position'); xlim([0 200])
ylabel('Trial #'); ylabel(cbar,'cm/s','FontSize',12,'Rotation',180)
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

sem = std(bnvel)/sqrt(sess.nlaps);
ciup = rmmissing(mean(bnvel) + sem*1.96);
cidn = rmmissing(mean(bnvel) - sem*1.96);

fhandle2 = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
plot(trackedges(1:end-1)*100,mean(bnvel),'k','LineWidth',2)
patch(100*[trackedges(1:length(cidn)),fliplr(trackedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
% plot(bnvel','Color',[.5 .5 .5])
xlabel('Position'); xlim([0 200])
ylabel('Average Velocity'); ylim([0 prctile(sess.velshft,99)*100])
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

end