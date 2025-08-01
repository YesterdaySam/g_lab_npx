function [binedges,binfr,ri,fhandle] = plot_frXrwdtime(root,unit,sess,tbnsz,wlen,plotflag)
%% Plots the avg binned firing rate by time to reward of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% tbnsz = size of position bins, default 0.05m = 5cm
% vthresh = threshold of behavioral velocity to throw out spikes, default 0.04 m/s
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
    tbnsz    = .25  %sec
    wlen     = 5    %sec On one side of the reward
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = -wlen:tbnsz:wlen;
spkinds = root.tsb(root.cl == unit);
% spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
spkts = sess.ts(spkinds);

spkmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpspks = spkts(spkts > sess.ts(sess.rwdind(i)) - wlen & spkts < sess.ts(sess.rwdind(i)) + wlen) - sess.ts(sess.rwdind(i));
    spkct = histcounts(tmpspks,binedges);
    % rtime   = sess.ts(sess.rwdind(i));
    % spkct   = histcounts(spkts, binedges + rtime);
    % runocc  = histcounts(sess.ts(sess.ind), binedges + rtime);
    % runocc  = runocc ./ sess.samprate;
    % bnoccs  = [bnoccs; runocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

bnoccs = ones(sess.nlaps,length(binedges)-1)*tbnsz;

% noOcc = bnoccs == 0;
% bnoccs(noOcc) = 0.0001; % Prevent div by 0 and NaNs later

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);

binfr = spksmooth ./ occsmooth;
rawfr = rmmissing(spkct ./ occct);

pOcc = sum(bnoccs) ./ sum(bnoccs,'all','omitnan');
uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
ri = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan');  % Reward Info in bits/spike

if plotflag
    sem = rmmissing(std(spkmap ./ bnoccs,'omitnan')/sqrt(sess.nlaps));
    ciup = rawfr + sem*1.96;
    cidn = rawfr - sem*1.96;
    % ciup = prctile(spkmap./bnoccs,95,1);
    % cidn = prctile(spkmap./bnoccs,5,1);

    binpos = binedges(1:end-1)+tbnsz/2;
    fhandle = figure; hold on
    plot(binpos,rawfr, 'k', 'LineWidth', 2)
    patch([binpos,fliplr(binpos)],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    plot(binpos,binfr, 'b')
    xlabel('Time (sec)'); ylabel('Firing Rate (spk/s)')

    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([0 10])
    elseif max(binfr,[],'all') < 20
        ylim([0 20])
    else
        ylim([0 max(rawfr)])
    end
    plot([0 0],[0 max(ylim)],'k--')

    % Plot velocity overlay
    ax = gca;
    yyaxis right
    if ~isfield(sess,'velXrwd')
        [~,bnvel] = get_velXrwd(sess,tbnsz,wlen,0);
    else
        bnvel = sess.velXrwd;
    end
    plot(binpos,bnvel,'r','LineWidth',2)
    ylim([0 prctile(sess.velshft,99)])
    ax.YAxis(2).Color = 'r';
    ylabel('Velocity (cm/s)')

    yyaxis left

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
end

end