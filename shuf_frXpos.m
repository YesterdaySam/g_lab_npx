function [binedges,binfr,fhandle] = shuf_frXpos(root,unit,sess,nShuf,dbnsz,plotflag)
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
% Created 5/28/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    nShuf = 100
    dbnsz = 0.05    %m
    plotflag = 1    %binary
end

% Only use valid trials by classifying error trials as non-running indices
for i = 1:length(sess.errTrials)
    sess.runInds(sess.lapstt(sess.errTrials(i)):sess.lapend(sess.errTrials(i))) = 0;
end

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods

spkct = histcounts(sess.pos(spkinds),binedges);
occct = histcounts(sess.pos(sess.runInds),binedges) / sess.samprate;

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);

binfr = spksmooth ./ occsmooth;
rawfr = rmmissing(spkct ./ occct);

for i = 1:nShuf
    [tmproot,tmpsess] = shiftTrain(root,sess);
    spkinds = tmproot.tsb(root.cl == unit);
    spkinds = spkinds(tmpsess.runInds(spkinds));   % Use only spikes in run periods

    tmpspkct(i,:) = histcounts(sess.pos(spkinds),binedges);
    tmpoccct = histcounts(sess.pos(tmpsess.runInds),binedges) / sess.samprate;

    tmpfr(i,:) = smoothdata(tmpspkct(i,:),'gaussian',5) ./ smoothdata(tmpoccct,'gaussian',5);
end

sem = rmmissing(std(tmpfr,'omitnan')/sqrt(nShuf));
ciup = mean(tmpfr) + sem*1.96;
cidn = mean(tmpfr) - sem*1.96;

if plotflag
    fhandle = figure; hold on
    plot([binedges(1:length(binfr))]*100,binfr, 'r')
    plot([binedges(1:length(binfr))]*100,mean(tmpfr), 'k')
    patch(100*[binedges(1:length(binfr)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    % plot([binedges(1:end-1)]*100,binfr, 'r')

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