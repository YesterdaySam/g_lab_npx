function [fhandle] = plot_trialHeatmap(root,unit,sess,dbnsz,vthresh,plotflag,normflag,smFactor)
%% Plots the spatially linearized spike raster for each trial of one unit
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   unit = cluster ID
%   sess = session struct from importBhvr
%   dbnsz = size of positional bins, default 0.03m = 2cm/s
%   vthresh = unused currently
%   plotflag = binary of whether to plot the output
%   normflag = whether to normalize heatmap values
%   smFactor = smoothing kernel for gaussian smoothdata on ratemap
%
% Outputs:
%   fhandle = handle to figure
%
% Created 7/10/24 LKW; Grienberger Lab; Brandeis University
% Updated to use get_frXpos core function 8/5/25 LKW
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %meters
    vthresh = 0.04  %meters/sec
    plotflag = 1    %binary
    normflag = 0    %binary
    smFactor = 5    %smoothing kernel size of rate map
end

% % Only use valid trials
% sess.lapstt = sess.lapstt(sess.valTrials);
% sess.lapend = sess.lapend(sess.valTrials);
% sess.nlaps  = length(sess.lapstt);
% 
% binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
% spkinds = root.tsb(root.cl == unit);
% % spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold
% spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
% spkmap = [];
% 
% for i = 1:sess.nlaps
%     tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
%     spkct   = histcounts(tmpspks, binedges);
%     bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges);
%     bnoccs  = bnocc / sess.samprate;
%     spkmap  = [spkmap; spkct ./ bnoccs];    % Normalize by time spent per bin for FR
% end

[binedges,spkmap] = get_frXpos(root,unit,sess,dbnsz,1.85,1,smFactor);

if plotflag
    fhandle = figure;
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
    colormap("parula")
    if normflag
        spkmap = normalize(spkmap,'range');
        imagesc(spkmap,[prctile(spkmap,1,'all'), prctile(spkmap,99,'all')]);
        xticks([1,length(binedges)-1]); xticklabels(binedges([1 end])*100);
        yticks(0:30:size(spkmap,1));
    else
        imagesc(spkmap,[prctile(spkmap,1,'all'), prctile(spkmap,99,'all')]);
        cbar = colorbar;
        ylabel(cbar,'FR (Hz)','FontSize',12,'Rotation',90)
        xticks(1:10:length(binedges)+1)
        xticklabels(binedges(1:10:end)*100)
    end
    clim([prctile(spkmap,1,'all'), prctile(spkmap,99,'all')]);
    xlabel('Position (cm)');
    ylabel('Trial #');
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end