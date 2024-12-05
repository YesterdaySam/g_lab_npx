function [binedges,binfr,fhandle] = plot_frXpos(root,unit,sess,dbnsz,vthresh,plotflag)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% dbnsz = size of position bins, default 0.05m = 5cm
% vthresh = threshold of behavioral velocity to throw out spikes, default 0.04 m/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% binedges = spatial bin edges
% binfr = spatial-binned firing rate
% fhandle = handle to figure
%
% Created 7/15/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %m
    vthresh = 0.04  %m/s; velocity threshold for spikes
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
% spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold 
spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods

% dspk = histcounts(sess.pos(spkinds),binedges);
% docc = histcounts(sess.pos,binedges)/sess.samprate;
% binfr = dspk ./ docc;

spkmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    spkct   = histcounts(tmpspks, binedges);
    lapInds = sess.ind(sess.lapstt(i):sess.lapend(i));
    tmpRun  = lapInds(sess.runInds(lapInds));   % Use only run periods for calculating occupancy
    % bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges) / sess.samprate;
    bnocc   = histcounts(sess.pos(tmpRun),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);

binfr = spksmooth ./ occsmooth;
rawfr = rmmissing(spkct ./ occct);

sem = rmmissing(std(spkmap ./ bnoccs,'omitnan')/sqrt(sess.nlaps));
ciup = rawfr + sem*1.96;
cidn = rawfr - sem*1.96;

if plotflag
    fhandle = figure; hold on
    plot([binedges(1:length(rawfr))]*100,rawfr, 'k')
    patch(100*[binedges(1:length(rawfr)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    plot([binedges(1:end-1)]*100,binfr, 'r')

    xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')

    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([0 10])
    elseif max(binfr,[],'all') < 20
        ylim([0 20])
    end

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
end

end