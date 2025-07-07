function [binedges,binfr,fhandle] = plot_frXpos(root,unit,sess,dbnsz,vFlag,plotflag)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% dbnsz = size of position bins, default 0.05m = 5cm
% vFlag = whether or not to remove spikes not coinciding with sess.runInds
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
    % vthresh = 0.04  %m/s; velocity threshold for spikes
    vFlag = 1
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);
sess.valTrials = 1:sess.nlaps;

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);

if vFlag
    % spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold
    spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
end

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
    xcoords = (binedges(1:end-1) + 0.5*dbnsz)*100;
    plot(xcoords,rawfr, 'k')
    patch([xcoords,fliplr(xcoords)],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(xcoords,binfr, 'b')
    
    xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')

    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([-1 10])
    elseif max(binfr,[],'all') < 20
        ylim([-1 20])
    end

    % Plot velocity overlay
    ax = gca;
    yyaxis right
    [~,bnvel] = plot_trialvel(sess,dbnsz,0);
    plot(xcoords,mean(bnvel),'r','LineWidth',2)
    ylim([0 prctile(sess.velshft,99)])
    ax.YAxis(2).Color = 'r';
    ylabel('Velocity (cm/s)')

    yyaxis left

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
end

end