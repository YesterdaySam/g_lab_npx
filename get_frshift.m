function [frShift,ps,frViol,crossLap,fhandle] = get_frshift(root,sess,coarseBin,fineBin,nTestBins,pThresh,plotflag)

arguments
    root
    sess
    coarseBin = 60    % time bin size (sec) over which to detect population-level z-score shifts
    fineBin   = 1     % time bin size( sec) over which to calculate avg lap FR
    nTestBins = 10    % number of bins on either side of a z-shift to use in a t-test
    pThresh   = 0.001 % t-test threshold very conservative
    plotflag  = 1
end

nCC = length(root.good);
nlap = length(sess.lapstt);
lapFRs = nan(nCC,nlap);
zrangeMat = [];

% Coarse for large Z-scored FR shifts across units
for i = 1:nCC
    [~,binedges,ccZ] = get_presence(root,root.good(i),sess,coarseBin);
    zrangeMat(i,:) = abs(diff(ccZ));
    binedges = binedges(1:end-2);
end

% Ignore high Z-score bins at edges of recording
hiZbins = find(median(zrangeMat) > 1);
if ~isempty(hiZbins)
    for i = 1:length(hiZbins)
        lappad(i,1) = sum(sess.ts(sess.lapstt) < binedges(hiZbins(i)));
        lappad(i,2) = sum(sess.ts(sess.lapstt) > binedges(hiZbins(i)));
    end

    laprm = logical(sum(lappad < nTestBins,2));
    hiZbins(laprm) = [];
    nZcross = length(hiZbins);
    ps =  nan(nCC,nZcross);

    % Fine grain FR per lap (could do this cleaner)
    for i = 1:nCC
        [ccFR,tmpedges] = get_presence(root,root.good(i),sess,fineBin);
        tmpedges = tmpedges(1:end-1);
        for j = 1:nlap
            lapFR = ccFR(tmpedges > sess.ts(sess.lapstt(j)) & tmpedges < sess.ts(sess.lapend(j)));
            lapFRs(i,j) = mean(lapFR,'omitnan');
        end
    end

    frShift = false(length(root.good),nZcross);

    % For each big delta in Z at pop level, t-test FR in nearby laps to find units in violation
    if nZcross > 0
        for i = 1:nZcross
            crossLap(i) = find(sess.ts(sess.lapstt) > binedges(hiZbins(i)),1);
            for j = 1:nCC
                [~,ps(j,i)] = ttest(lapFRs(j,crossLap(i)-nTestBins:crossLap(i)-1),lapFRs(j,crossLap(i):crossLap(i)+nTestBins-1));
            end
            frShift(ps(:,i) < pThresh,i) = true;
            frViol(i) = sum(frShift(:,i)) / length(root.good);
        end
    else
        frViol = 0;
        crossLap = [];
    end

else
    ps = nan(nCC,1);
    frShift = false(nCC,1);
    frViol = 0;
    crossLap = [];
end

if plotflag
    uZxlap = mean(zrangeMat,'omitnan');
    normUZ = uZxlap ./ max(uZxlap) * nCC;
    fhandle = figure; hold on;
    imagesc(binedges,1:length(root.good),zrangeMat)
    set(gca, 'YDir', 'normal')
    colormap('hot')
    clim([0 5])
    ylabel("Unit #"); ylim([0 length(root.good)])
    xlabel('Time (s)'); xlim([0 binedges(end)])
    plot(binedges, normUZ,'w');
    plot(binedges(hiZbins),normUZ(hiZbins),'w*')
end

end
