function [binedges,binfr,fhandle] = plot_AligntoRand(root,unit,sess,dbnsz,vthresh,plotflag)
%% Plots the avg binned firing rate around the random cue
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
% FG 12/12
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %m
    vthresh = 0.04  %m/s; velocity threshold for spikes
    plotflag = 1    %binary
end

% probably fine since only looking around rand
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = -0.8:dbnsz:0.8;    % 100 cm (80 on each side of rand cue)
spkinds = root.tsb(root.cl == unit);
% spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold 
spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods

% dspk = histcounts(sess.pos(spkinds),binedges);
% docc = histcounts(sess.pos,binedges)/sess.samprate;
% binfr = dspk ./ docc;

spkmap = [];
bnoccs = [];
beltLength = 1.8;
for i = 2:length(sess.randCueInd)-1
    
    if sess.pos(sess.randCueInd(i))-0.8 < 0
        beforeCue = beltLength + (sess.pos(sess.randCueInd(i)) - 0.8);
    else
        beforeCue = sess.pos(sess.randCueInd(i))-0.8;
    end
    
    if sess.pos(sess.randCueInd(i))+0.8 > beltLength
        afterCue = abs(beltLength - sess.pos(sess.randCueInd(i))+0.8);
    else
        afterCue = sess.pos(sess.randCueInd(i))+0.8;
    end
    
    %something wrong
    try
        initInd = find(sess.pos(1:sess.randCueInd(i)) <= beforeCue+0.01);
        initInd = initInd(end);
    catch
    end
    lastInd = find(sess.pos(sess.randCueInd(i):end) >= afterCue);
    lastInd = lastInd(1) + sess.randCueInd(i);  

    tmpspks = sess.pos(spkinds(spkinds > initInd & spkinds < lastInd)) - sess.pos(sess.randCueInd(i));
    spkct   = histcounts(tmpspks, binedges);

    lapInds = sess.ind(initInd:lastInd);  

    tmpRun  = lapInds(sess.runInds(lapInds));   % Use only run periods for calculating occupancy
    % bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges) / sess.samprate;
    bnocc   = histcounts(sess.pos(tmpRun)-sess.pos(sess.randCueInd(i)),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);

%needs NaN to be 0
binfr = spksmooth ./ occsmooth;
binfr(isnan(binfr)) = 0;
rawfr = spkct ./ occct;
rawfr(isnan(rawfr)) = 0;

sem = std(spkmap ./ bnoccs,'omitnan')/sqrt(sess.nlaps);
sem(isnan(sem)) = 0;
ciup = rawfr + sem*1.96; %weird for some laps big SEM?
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
    xline(0)
end

end