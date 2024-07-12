function [fhandle1, fhandle2] = plot_lickraster_rwd(sess,tsrange,bnsz)
%% plot a peri-reward lick raster given a session struct
% Inputs
% session = struct from importBhvr.m
% tsrange = double of padding (in seconds) on either side of reward
% bnsz    = double of bin size (in seconds) for lick rate graph (e.g. 0.1)
% Outputs
% fhandle = handle to figure

arguments
    sess
    tsrange = 5     % seconds of padding
    bnsz = 0.1      % seconds
end

nrwd    = length(sess.rwdind);
tsbuff  = tsrange * sess.samprate;
lckmap  = [];

for i = 1:nrwd
    tmprwd = sess.rwdind(i);
    tmpstt = tmprwd - tsbuff;
    tmpend = tmprwd + tsbuff;
    tmplck = find(sess.lckind > tmpstt & sess.lckind < tmpend);
    lckt   = sess.ts(sess.lckind(tmplck))';
    lckmap = [lckmap; lckt - sess.ts(tmprwd), i*ones(numel(tmplck),1)];
    pslick(i,:) = histcounts(sess.lckind(tmplck),tmpstt:bnsz*sess.samprate:tmpend);
end

fhandle1 = figure;      % Lick Raster
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
% plot(lckmap(:,1),lckmap(:,2),'|','Color',[.65 .65 1])
plot(lckmap(:,1),lckmap(:,2),'k|')
xlabel('Time to reward release')
ylabel('Trial #')
set(gca,'FontSize',12,'FontName','Arial')
title([strtok(sess.name,'_'), ' ', sess.name(end-9:end-8)])

fhandle2 = figure;      % Smoothed average lickrate
set(gcf,'units','normalized','position',[0.4 0.1 0.3 0.2])
plot(-tsrange+bnsz:bnsz:tsrange, mean(pslick)/bnsz)
xlabel('Time to reward release')
ylabel('Avg. Licks/sec')
title([strtok(sess.name,'_'), ' ', sess.name(end-9:end-8)])
set(gca,'FontSize',12,'FontName','Arial')

end