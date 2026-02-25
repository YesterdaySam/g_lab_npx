function [fhandle] = plot_burstSpkRaster(root,unit,sess,vFlag,plotflag)
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
    vFlag = 1  %meters/s; velocity threshold
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

spkinds = root.tsb(root.cl == unit);
bstinds = get_bursts(root,sess,unit);
if vFlag
    spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
    bstinds = bstinds(sess.runInds(bstinds));
end

spkpos = [];
bstpos = [];
rwdpos = [];

for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    tmpbsts = sess.pos(bstinds(bstinds > sess.lapstt(i) & bstinds < sess.lapend(i)));
    tmprwd  = sess.pos(sess.rwdind(sess.rwdind > sess.lapstt(i) & sess.rwdind < sess.lapend(i)));

    spkpos  = [spkpos; tmpspks, i*ones(length(tmpspks),1)];
    bstpos  = [bstpos; tmpbsts, i*ones(length(tmpbsts),1)];

    if isempty(tmprwd)
        rwdpos = [rwdpos; NaN, i];
    else
        rwdpos  = [rwdpos; tmprwd, i];
    end
end

% spkpos = sess.pos(root.tsb(root.cl == unit));

if plotflag
    fhandle = figure;      % Positional Lick Raster
    hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.21 0.36])
    plot(spkpos(:,1)*100,spkpos(:,2),'k|')
    plot(bstpos(:,1)*100,bstpos(:,2),'ro')
    plot(rwdpos(:,1)*100,rwdpos(:,2),'b*')
    xlabel('Position (cm)'); xlim([0 100*max(sess.pos(sess.lapstt(1):sess.lapend(1)))])
    ylabel('Trial #'); ylim([0 sess.nlaps])
    set(gca,'FontSize',12,'FontName','Arial')
end

end