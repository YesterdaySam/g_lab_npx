function [fhandle1, fhandle2] = plot_lickraster_rwd(session,tsrange,bnsz)
%% plot a peri-reward lick raster given a session struct
% Inputs
% session = struct from importBhvr.m
% tsrange = double of padding (in seconds) on either side of reward
% bnsz    = double of bin size (in seconds) for lick rate graph (e.g. 0.1)
% Outputs
% fhandle = handle to figure

nrwd    = length(session.rwdind);
tsbuff  = tsrange * session.samprate;
lckmap  = [];

for i = 1:nrwd
    tmprwd = session.rwdind(i);
    tmpstt = tmprwd - tsbuff;
    tmpend = tmprwd + tsbuff;
    tmplck = find(session.lckind > tmpstt & session.lckind < tmpend);
    lckt   = session.ts(session.lckind(tmplck))';
    lckmap = [lckmap; lckt - session.ts(tmprwd), i*ones(numel(tmplck),1)];
    pslick(i,:) = histcounts(session.lckind(tmplck),tmpstt:bnsz*session.samprate:tmpend);
end

fhandle1 = figure;      % Lick Raster
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
% plot(lckmap(:,1),lckmap(:,2),'|','Color',[.65 .65 1])
plot(lckmap(:,1),lckmap(:,2),'k|')
xlabel('Time to reward release')
ylabel('Trial #')
set(gca,'FontSize',12,'FontName','Arial')
title([strtok(session.name,'_'), ' ', session.name(end-9:end-8)])

fhandle2 = figure;      % Smoothed average lickrate
set(gcf,'units','normalized','position',[0.4 0.1 0.3 0.2])
plot(-tsrange+bnsz:bnsz:tsrange, mean(pslick)/bnsz)
xlabel('Time to reward release')
ylabel('Avg. Licks/sec')
title([strtok(session.name,'_'), ' ', session.name(end-9:end-8)])
set(gca,'FontSize',12,'FontName','Arial')

end