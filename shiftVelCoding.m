function [shiftRs,pkShiftTime,pkCorr,fhandle] = shiftVelCoding(root,unit,sess,winlen,vbnsz,corrType,plotflag)

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    winlen = 1      % sec
    vbnsz = 2       %cm/s
    corrType = 'Pearson'
    plotflag = 0    %binary
end

shifts = -winlen*sess.samprate:10:winlen*sess.samprate;
nShifts = length(shifts);

binedges = 0:vbnsz:max(sess.velshft);

spktsb = root.tsb(root.cl == unit);

shiftRs = zeros(nShifts,1);

vocc = histcounts(sess.velshft,binedges)/sess.samprate;
lowOcc = vocc < 1; % Ignore bins with <1sec occupancy
lowOcc(1:2) = logical([1 1]); % Exclude periods <4cm/s velocity
usebins = logical([0 1-lowOcc(1:end)]);  %Drop first 2 bins (4cm/s) of standing periods and low occ bins
binctrs = binedges(usebins) - vbnsz/2;

for i = 1:nShifts
    shiftSpks = spktsb + shifts(i); % Indices
    cleanInds = shiftSpks < 0 | shiftSpks > sess.ind(end);
    shiftSpks(cleanInds) = []; % Remove unreachable spikes
    shiftSpkVel = sess.velshft(shiftSpks);

    vspk = histcounts(shiftSpkVel,binedges);
    binfr = vspk ./ vocc;
    binfr(lowOcc) = [];

    tmpCorr = corr(binctrs',binfr','Type',corrType);
    shiftRs(i) = tmpCorr;
end

[pkCorr,pkloc] = max(shiftRs);
[mnCorr,mnloc] = min(shiftRs);
if abs(mnCorr) > pkCorr
    pkCorr = mnCorr;
    pkloc = mnloc;
end
pkShiftTime = shifts(pkloc)/sess.samprate;


[ccorr,clags] = xcorr(sess.velshft,histcounts(sess.ts(spktsb),sess.ts),winlen*sess.samprate);
pkShiftTime = clags(ccorr == max(ccorr));


if plotflag

    fhandle = figure;
    set(gcf,'units','normalized','position',[0.4 0.15 0.2 0.7])

    subplot(2,1,1); hold on
    plot([0 0],[-1 1],'k--',"HandleVisibility",'off')
    plot(shifts/sess.samprate,shiftRs,'k',"HandleVisibility",'off')
    plot(shifts(pkloc)/sess.samprate,pkCorr,'r*')
    xlabel('Lag (s)'); ylabel(['R (' corrType ')'])
    legend({'Best lag'},'Location','best')
    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')

    subplot(2,1,2); hold on
    lag0Spks = sess.velshft(spktsb);
    bestSpks = spktsb + shifts(pkloc);
    cleanInds = bestSpks < 0 | bestSpks > sess.ind(end);
    bestSpks(cleanInds) = []; % Remove unreachable spikes
    shiftSpkVel = sess.velshft(bestSpks);
    vspkbest = histcounts(shiftSpkVel,binedges);
    vspk0 = histcounts(lag0Spks,binedges);
    binfrbest = vspkbest ./ vocc;
    binfrbest(lowOcc) = [];
    binfr0 = vspk0 ./ vocc;
    binfr0(lowOcc) = [];
    corrBest = corr(binctrs',binfrbest');
    corr0 = corr(binctrs',binfr0');
    
    scatter(binctrs,binfrbest,'k','filled')
    % scatter(binctrs,binfr0,repmat([.5 .5 .5],[length(binfr0),1]),'filled')
    scatter(binctrs,binfr0,'k')
    xlabel('Velocity (cm/s)'); ylabel('Firing Rate (spk/s)')
    set(gca,'FontSize',12,'FontName','Arial')
    legend({'Best lag','Lag 0'},'Location','sw')
    ylims = ylim;
    ylim([0 ylims(2)+ylims(2)*0.3])
    ylims = ylim;
    xlims = xlim;
    xlim([0 xlims(2)])
    xlims = xlim;

    if corrBest > 0
        text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims), ['R Lag = ' num2str(corrBest, 5)], 'FontSize', 12)
        text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['Lag (ms) = ' num2str(1000*shifts(pkloc)/sess.samprate, 3)], 'FontSize', 12)
        text(xlims(2) - .9*diff(xlims), ylims(2)-.2*diff(ylims), ['R 0 = ' num2str(corr0, 5)], 'FontSize', 12)
    else
        text(xlims(2) - .5*diff(xlims), ylims(2)-.1*diff(ylims), ['R Lag = ' num2str(corrBest, 5)], 'FontSize', 12)
        text(xlims(2) - .5*diff(xlims), ylims(2)-.15*diff(ylims), ['Lag (ms) = ' num2str(1000*shifts(pkloc)/sess.samprate, 3)], 'FontSize', 12)
        text(xlims(2) - .5*diff(xlims), ylims(2)-.2*diff(ylims), [ 'R 0 = ' num2str(corr0, 5)], 'FontSize', 12)
    end

end

end