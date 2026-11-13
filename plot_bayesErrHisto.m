function [errLocs,errMean, fhandle] = plot_bayesErrHisto(binedges,decode1,decode2,decode3,plotflag)

arguments
    binedges
    decode1
    decode2
    decode3
    plotflag = 1
end

[loc1,pErr1] = get_decodeErrLoc(decode1.dErr,binedges);
[loc2,pErr2] = get_decodeErrLoc(decode2.dErr,binedges);
[loc3,pErr3] = get_decodeErrLoc(decode3.dErr,binedges);

binpos = binedges(2:end) - 0.5*diff(binedges(1:2));

errLocs = [loc1, loc2, loc3];
errMean = [mean(abs(decode1.dErr)), mean(abs(decode2.dErr)), mean(abs(decode3.dErr))];

if plotflag
    fhandle = figure; hold on;
    plot(binpos, pErr1, 'k')
    plot(binpos, pErr2, 'r')
    plot(binpos, pErr3, 'c')
    plot(loc1, max(pErr1)+0.005, 'kv')
    plot(loc2, max(pErr2)+0.005, 'rv')
    plot(loc3, max(pErr3)+0.005, 'cv')

    xlabel('Abs. error (m)')
    ylabel('Probability')
    set(gca,'FontSize',16,'FontName','Arial')

    legend('F by F', 'N by N', 'N by F')
end
end