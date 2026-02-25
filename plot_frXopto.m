function [bnctrs,pulseMat,fhandle] = plot_frXopto(root,unit,sess,bnsz,bnrg,plotflag)
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
% Created 9/23/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    bnsz = 0.002    %sec
    bnrg = 0.5      %sec on either side of the pulse start
    plotflag = 1    %binary
end

nPulse = length(sess.optoind);
binedges = -bnrg:bnsz:bnrg;
nBins = length(binedges)-1;
bnctrs = -bnrg+0.5*bnsz:bnsz:bnrg;

spkinds = root.tsb(root.cl == unit);
spkts   = sess.ts(spkinds);

pulseMat = zeros(nPulse,nBins);

for i = 1:nPulse
    tmpT = sess.ts(sess.optoind(i));
    tmpedges = tmpT-bnrg:bnsz:tmpT+bnrg;
    pulseMat(i,:) = histcounts(spkts,tmpedges);
end

if plotflag
    fhandle = figure; hold on
    b1 = bar(bnctrs*1000,sum(pulseMat),'FaceColor',[0.5 0.5 1]);
    plot(bnctrs*1000,smooth(sum(pulseMat),5),'r');
    lims = ylim;
    plot([0 0], lims, 'k--')
    xlabel('Time after opto pulse (ms)'); ylabel('Spike count')
    % histogram(pulseMat,-bnrg:bnsz:bnrg,'Normalization','probability')
    % plot([binedges(1:end-1)]*100,mean(binfr,1,'omitnan'), 'k')
    % patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    % 
    % xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')
    % 
    % if max(mean(binfr,1,'omitnan'),[],'all') < 10
    %     ylim([0 10])
    % elseif max(binfr,[],'all') < 20
    %     ylim([0 20])
    % end

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
end

end