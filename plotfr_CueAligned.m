function [binedges,binfr,fhandle] = plotfr_CueAligned(root,unit,sess,cueFlag,dbnsz,vFlag,plotflag)
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
    cueFlag         %1 for fixed, 0 for rand
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

%find bin of cue
if cueFlag %align to fixed cue
    posCue = sess.pos(sess.fixedCueInd(1));
    binCue = find(binedges(:)>posCue);
    binCue = binCue(1)-1; %cue is in this bin
    %circshift for binCue to be at center
    center = floor(length(binedges)/2)+1;
    moveby = center - binCue; 

    spkmap = circshift(spkmap, [0 moveby]);
    bnoccs = circshift(bnoccs, [0 moveby]);
else %align to rand cue
    newbnoccs = [];
    newspkmap = [];
    
    if cueFlag
        maxNum = sess.nlaps-1;
    else
        maxNum = length(sess.randCueInd);
        if sess.nlaps-1 < length(sess.randCueInd)
            maxNum = sess.nlaps-1;
        end
        %maxNum = sess.nlaps/2;
    end

    for j = 2:maxNum 
        currSpkMap = spkmap(j,:);
        currBnoccs = bnoccs(j,:);

        posCue = sess.pos(sess.randCueInd(j));
        binCue = find(binedges(:)>posCue);
        if isempty(binCue)
            continue;
        end
        binCue = binCue(1)-1; %cue is in this bin
        %circshift for binCue to be at center
        % center = floor(length(binedges)/2)+1;
        % moveby = center - binCue; 

        % currSpkMap = circshift(currSpkMap, moveby);
        % currBnoccs = circshift(currBnoccs, moveby);
        binStart = binCue - 19;
        binEnd = binCue + 19;

        if binStart < 0 
            prevLapSpkMap = spkmap(j-1,:);
            prevLapBnoccs = bnoccs(j-1,:);
            tempNewSpkMap = prevLapSpkMap(37+binStart:end);
            tempNewBnoccs = prevLapBnoccs(37+binStart:end);

            leftOver = 19 + binStart;
            tempNewSpkMap = [tempNewSpkMap currSpkMap(1:leftOver-1)];
            tempNewBnoccs = [tempNewBnoccs currBnoccs(1:leftOver-1)];
        else
            tempNewSpkMap = currSpkMap(binStart+1:binCue);
            tempNewBnoccs = currBnoccs(binStart+1:binCue);
        end

        if binEnd > 37
            tempNewSpkMap = [tempNewSpkMap currSpkMap(binCue+1:37)];
            tempNewBnoccs = [tempNewBnoccs currBnoccs(binCue+1:37)];

            leftOver = binEnd - 37;
            nextLapSpkMap = spkmap(j+1,:);
            nextLapBnoccs = bnoccs(j+1,:);

            tempNewSpkMap = [tempNewSpkMap nextLapSpkMap(1:leftOver-1)];
            tempNewBnoccs = [tempNewBnoccs nextLapBnoccs(1:leftOver-1)];
        else 
            tempNewSpkMap = [tempNewSpkMap currSpkMap(binCue+1:binEnd-1)];
            tempNewBnoccs = [tempNewBnoccs currBnoccs(binCue+1:binEnd-1)];
        end

        newbnoccs  = [newbnoccs; tempNewBnoccs];              % Save bin occupancy
        newspkmap  = [newspkmap ; tempNewSpkMap]; 
    end
    spkmap = newspkmap;
    bnoccs = newbnoccs;
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
    % 
    % % Plot velocity overlay
    % ax = gca;
    % yyaxis right
    % if ~isfield(sess,'sess.velXpos')
    %     [~,bnvel] = plot_trialvel(sess,dbnsz,0);
    % else
    %     bnvel = sess.velXpos;
    % end
    % velmean = mean(bnvel);
    % velsem = std(bnvel)./sqrt(sess.nlaps);

    % plot(xcoords,velmean,'r','LineWidth',2)
    % patch([xcoords,fliplr(xcoords)],[velmean-velsem*1.96,fliplr(velmean+velsem*1.96)],'r','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    % ylim([0 prctile(sess.velshft,99)])
    % ax.YAxis(2).Color = 'r';
    % ylabel('Velocity (cm/s)')

    xcoords = xcoords - xcoords(19);

    plot(xcoords, rawfr, 'k-')
    patch([xcoords,fliplr(xcoords)],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
    plot(xcoords, binfr, 'b-')
    
    xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')

    if cueFlag
        title(['Fixed Cue Aligned ' 'Unit ' num2str(unit)])
    else
        title(['Rand Cue Aligned ' 'Unit ' num2str(unit)])
    end

    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([-1 10])
    elseif max(binfr,[],'all') < 20
        ylim([-1 20])
    end
    xline(0)

    set(gca,'FontSize',12,'FontName','Arial')
end

end