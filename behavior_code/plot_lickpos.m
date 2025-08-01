function [binedges, lckmap, bnlck, fhandle, fhandle2] = plot_lickpos(sess,bnsz,plotflag)
%% Create linearized lick posiition (binned by space)
% Inputs
%   sess = struct from importBhvr.m
%   bnsz    = double in meters (m)
%
% Outputs
%   fhandle = handle to figure
%
% Created 5/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    bnsz = 0.03      % in m
    plotflag = 1
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
nlaps       = length(sess.lapstt);

binedges = 0:bnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial

nbins           = size(binedges,2);

% Find licks
lckmap  = [];
bnlck = ones(nlaps,nbins-1);

for i = 1:nlaps
    % Find index of licks for this lap
    tmplck = sess.lckind(find(sess.lckind > sess.lapstt(i) & sess.lckind < sess.lapend(i)));
    
    % Create Nx2 of [lick times, trial]
    lckmap = [lckmap; sess.pos(tmplck), i*ones(numel(tmplck),1)];
    
    bnocc = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges);
    bnlck(i,:) = histcounts(sess.pos(tmplck),binedges);
    bnlck(i,:) = bnlck(i,:) ./ (bnocc / sess.samprate);  %Normalize to time in each bin
end

% Find rewards after valid lap starts
rwdmap = [];
for i = 1:nlaps
    tmprwd = sess.pos(sess.rwdind(find(sess.rwdind > sess.lapstt(i),1)));
    rwdmap = [rwdmap; tmprwd, i];
end

if plotflag
    fhandle = figure;      % Positional Lick Raster
    hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    plot(lckmap(:,1)*100,lckmap(:,2),'k|')
    plot(rwdmap(:,1)*100,rwdmap(:,2),'r*')
    xlabel('Position (cm)')
    ylabel('Trial #')
    set(gca,'FontSize',12,'FontName','Arial')

    sem = std(bnlck,'omitnan')/sqrt(nlaps);
    ciup = rmmissing(mean(bnlck,1,'omitnan') + sem*1.96);
    cidn = rmmissing(mean(bnlck,1,'omitnan') - sem*1.96);

    fhandle2 = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
    plot(binedges(1:end-1)*100,mean(bnlck,1,'omitnan'),'b','LineWidth',2)
    patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'b','FaceAlpha',0.5,'EdgeColor','none')
    xlabel('Position (cm)'); xlim([0 200])
    ylabel('Average Licks/s');
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end