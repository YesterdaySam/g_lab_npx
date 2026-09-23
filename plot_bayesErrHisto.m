function [errLocs,errMean, fhandle] = plot_bayesErrHisto(binedges,decode1,decode2,decode3)

errDist1 = histcounts(abs(decode1.dErr),binedges);
errDist2 = histcounts(abs(decode2.dErr),binedges);
errDist3 = histcounts(abs(decode3.dErr),binedges);

binpos = binedges(2:end) - 0.5*diff(binedges(1:2));

[~,loc1] = max(errDist1);
[~,loc2] = max(errDist2);
[~,loc3] = max(errDist3);

pErr1 = errDist1 ./ sum(errDist1);
pErr2 = errDist2 ./ sum(errDist2);
pErr3 = errDist3 ./ sum(errDist3);
errLocs = binpos([loc1, loc2, loc3]);
errMean = [mean(abs(decode1.dErr)), mean(abs(decode2.dErr)), mean(abs(decode3.dErr))];

fhandle = figure; hold on;
plot(binpos, pErr1, 'k')
plot(binpos, pErr2, 'r')
plot(binpos, pErr3, 'c')
plot(binpos(loc1), pErr1(loc1)+0.005, 'kv')
plot(binpos(loc2), pErr2(loc2)+0.005, 'rv')
plot(binpos(loc3), pErr3(loc3)+0.005, 'cv')

xlabel('Abs. error (m)')
ylabel('Probability')
set(gca,'FontSize',16,'FontName','Arial')

legend('F by F', 'N by N', 'N by F')

end