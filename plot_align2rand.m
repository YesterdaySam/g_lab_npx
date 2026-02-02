function [binedges,binfr,fhandle,rastfig] = plot_align2rand(root,unit,sess,qflag,dbnsz,plotflag)
%% Plots the avg binned firing rate around the random cue
% For speed, pre-calculate random cue starts with get_QTwin.m function
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   unit = cluster ID
%   sess = session struct from importBhvr
%   dbnsz = size of position bins, default 0.05m = 5cm
%   plotflag = binary of whether to plot the output
%
% Outputs:
%   binedges = spatial bin edges
%   binfr = spatial-binned firing rate
%   fhandle = handle to figure
%
% Created 1/27/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    qflag = 0       % 0 = rand cue; 1 = fixed cue
    dbnsz = 0.05    %m
    plotflag = 1    %binary
end

beltLength = 1.85;  % Hard coded for now
wlen = beltLength / 2; % meters
binedges = -wlen:dbnsz:wlen;    % 100 cm (80 on each side of rand cue)
spkinds = root.tsb(root.cl == unit);
spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods

spkmap = [];
bnoccs = [];
rastMt = [];

if qflag == 0
    qs = sess.randCueInd;
    stts = sess.rqStt;
    ends = sess.rqEnd;
else 
    qs = sess.fixedCueInd;
    stts = sess.fqStt;
    ends = sess.fqEnd;
end

for i = 2:length(qs)-1
    tmpCue = qs(i);    % Set to sess.fixedCueInd to validate against plot_frXpos.m trace
    tmplap = find(sess.lapstt < tmpCue,1,'last');

    % Find spikes, align to rand cue and bin
    tmpspks = spkinds(spkinds > stts(i) & spkinds < ends(i));
    pstLapSpks = tmpspks > sess.lapstt(tmplap+1);
    preLapSpks = tmpspks < sess.lapend(tmplap-1);
    tmpspkpos = sess.pos(tmpspks) - sess.pos(tmpCue);
    tmpspkpos(pstLapSpks) = tmpspkpos(pstLapSpks) + beltLength;
    tmpspkpos(preLapSpks) = tmpspkpos(preLapSpks) - beltLength;
    
    % Calculate occupancy over run periods and align to rand cue
    lapInds = sess.ind(stts(i):ends(i));  
    tmpRun  = lapInds(sess.runInds(lapInds));
    pstLapRun = tmpRun > sess.lapstt(tmplap+1);
    preLapRun = tmpRun < sess.lapend(tmplap-1);
    tmpPos = sess.pos(tmpRun) - sess.pos(tmpCue);
    tmpPos(pstLapRun) = tmpPos(pstLapRun) + beltLength;
    tmpPos(preLapRun) = tmpPos(preLapRun) - beltLength;

    % Store lap info
    bnoccs = [bnoccs; histcounts(tmpPos,binedges) / sess.samprate];
    spkmap = [spkmap; histcounts(tmpspkpos,binedges)];
    rastMt = [rastMt; tmpspkpos tmplap*ones(size(tmpspkpos))];
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

    if qflag == 0
        xlabel('Distance to Rand Cue (cm)'); ylabel('Firing Rate (spk/s)')
    else
        xlabel('Distance to Fixed Cue (cm)'); ylabel('Firing Rate (spk/s)')
    end

    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([0 10])
    elseif max(binfr,[],'all') < 20
        ylim([0 20])
    end

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
    
    % Raster
    rastfig = figure; plot(rastMt(:,1)*100,rastMt(:,2),'k|')
    if qflag == 0
        xlabel('Distance to Rand Cue (cm)'); ylabel('Trial')
    else
        xlabel('Distance to Fixed Cue (cm)'); ylabel('Trial')
    end
    set(gca,'FontSize',12,'FontName','Arial')
end

end