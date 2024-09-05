function [fhandle] = plot_trialHeatmap(root,unit,sess,dbnsz,vthresh,plotflag)
%% Plots the spatially linearized spike raster for each trial of one unit
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   unit = cluster ID
%   sess = session struct from importBhvr
%   dbnsz = size of positional bins, default 0.03m = 2cm/s
%   plotflag = binary of whether to plot the output
%
% Outputs:
%   fhandle = handle to figure
%
% Created 7/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.03    %meters
    vthresh = 0.02  %meters/sec
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold
spkmap = [];

for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    spkct   = histcounts(tmpspks, binedges);
    bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges);
    bnoccs  = bnocc / sess.samprate;
    spkmap  = [spkmap; spkct ./ bnoccs];    % Normalize by time spent per bin for FR
end

if plotflag
    fhandle = figure;
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    imagesc(spkmap,[prctile(spkmap,1,'all'), prctile(spkmap,98,'all')]);
    colormap("parula")
    cbar = colorbar; clim([0 inf]);
    xticks(1:10:length(binedges))
    xticklabels(binedges(1:10:end)*100)
    xlabel('Position (cm)');
    ylabel('Trial #'); ylabel(cbar,'FR (Hz)','FontSize',12,'Rotation',90)
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end