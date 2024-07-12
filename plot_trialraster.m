function [fhandle] = plot_trialraster(root,unit,sess,dbnsz,plotflag)
%% Plots the binned firing rate by velocity of a unit
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
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

% binedges = 0:dbnsz:max(sess.pos);
spkinds = root.tsb(root.cl == unit);

spkpos = [];
rwdpos = [];

for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    tmprwd  = sess.pos(sess.rwdind(sess.rwdind > sess.lapstt(i) & sess.rwdind < sess.lapend(i)));
    spkpos  = [spkpos; tmpspks, i*ones(length(tmpspks),1)];
    rwdpos  = [rwdpos; tmprwd, i];
end

% spkpos = sess.pos(root.tsb(root.cl == unit));

if plotflag
    fhandle = figure;      % Positional Lick Raster
    hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    plot(spkpos(:,1)*100,spkpos(:,2),'k.')
    plot(rwdpos(:,1)*100,rwdpos(:,2),'b*')
    xlabel('Position (cm)')
    ylabel('Trial #')
    set(gca,'FontSize',12,'FontName','Arial')
end

end