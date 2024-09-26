function [bnctrs,pulseMat,fhandle] = plot_frXopto_shift(root,ctl,unit,sess,bnsz,bnrg,plotflag)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% bnsz = size of histogram bins, default 0.02m/s = 2cm/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% binedges = spatial bin edges
% binfr = spatial-binned firing rate
% fhandle = handle to figure
%
% Created 9/24/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    ctl
    unit {double}   %Cluster ID
    sess            %session struct
    bnsz = 0.002    %sec
    bnrg = 0.1     %sec on either side of the pulse start
    plotflag = 1    %binary
end

nPulse = length(ctl);
binedges = -bnrg:bnsz:bnrg;
nBins = length(binedges)-1;

spkinds = root.tsb(root.cl == unit);
spkts   = sess.ts(spkinds);

pulseMat = zeros(nPulse,nBins);

for i = 1:nPulse
    tmpT = ctl(i);
    tmpedges = tmpT-bnrg:bnsz:tmpT+bnrg;
    pulseMat(i,:) = histcounts(spkts,tmpedges);
end

if plotflag
    fhandle = figure; hold on
    set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.2])
    bnctrs = -bnrg+0.5*bnsz:bnsz:bnrg;
    b1 = bar(bnctrs,sum(pulseMat),'FaceColor',[0.5 0.5 1]);
    plot(bnctrs,smooth(sum(pulseMat),5),'r');
    lims = ylim;
    plot([0 0], lims, 'k--')
    xlabel('Time after opto pulse (sec)'); ylabel('Spike count')
    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
else
    bnctrs = -bnrg+0.5*bnsz:bnsz:bnrg;
end

end