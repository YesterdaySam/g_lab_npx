function [fhandle] = plot_trialvel(session, bnsz)
%% Create linearized velocity (binned by space)
% Inputs
% session = struct from importBhvr.m
% bnsz    = double in cm e.g. 1
% Outputs
% fhandle = handle to figure

% vel             = session.vel - min(session.vel);
vel             = session.velshft;
bnsz            = bnsz/100;    % translate to m
tracklen        = (max(session.pos) - min(session.pos)); % m
trackedges      = 0:bnsz:tracklen;
nbins           = size(trackedges,2);
clear bnvel

for i = 1:session.nlaps-1
    bnocc = histcounts(session.pos(session.lapstt(i):session.lapend(i+1)), trackedges);
    for j = 1:nbins-1
        tmpind = session.lapstt(i) + find(session.pos(session.lapstt(i):session.lapend(i+1)) > trackedges(j) & session.pos(session.lapstt(i):session.lapend(i+1)) < trackedges(j+1));
        bnvel(i,j) = mean(vel(tmpind));
    end
end

fhandle = figure;
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
imagesc(bnvel,[prctile(bnvel,1,'all'), prctile(bnvel,99,'all')]); 
% imagesc(bnvel, [0 max(bnvel,[],'all')])
colormap("sky")
colorbar; clim([0 inf]);
xlabel('Position')
ylabel('Trial #')
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
title([strtok(session.name,'_'), ' ', session.name(end-9:end-8)])

end