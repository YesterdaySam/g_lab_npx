function [binedges,binfr,fhandle] = plot_frXpos(root,unit,sess,dbnsz,vthresh,plotflag)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% vbnsz = size of velocity bins, default 0.02m/s = 2cm/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% binedges = velocity bin edges
% binfr = velocity-binned firing rate
% mdlparams = R and p values of the correlation coefficient, slope and
%   y-intercept of a linear model fit between velocity and firing rate
% fhandle = handle to figure
%
% Created 7/15/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %m
    vthresh = 0.04  %m/s; velocity threshold for spikes
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold

% dspk = histcounts(sess.pos(spkinds),binedges);
% docc = histcounts(sess.pos,binedges)/sess.samprate;
% binfr = dspk ./ docc;

spkmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    spkct   = histcounts(tmpspks, binedges);
    bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

binfr = spkmap ./ bnoccs;

sem = std(binfr,'omitnan')/sqrt(sess.nlaps);
ciup = rmmissing(mean(binfr,1,'omitnan') + sem*1.96);
cidn = rmmissing(mean(binfr,1,'omitnan') - sem*1.96);

if plotflag
    fhandle = figure; hold on
    plot([binedges(1:end-1)]*100,mean(binfr,1,'omitnan'), 'k')
    patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')

    xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')

    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([0 10])
    elseif max(binfr,[],'all') < 20
        ylim([0 20])
    end

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
end

%% Spatial information
% SI = sum over all bins ( normalized occupancy .* FRmap .* log2 (FRmap ./ uFR) ) ./ uFR 
% FRmap is the number of spikes in each spatial bin
% normalized occupancy is the probability of the mouse in each spatial bin
% uFR is the overall average firing rate of the unit

uFR = sum(spkmap,'all','omitnan') / sum(bnoccs,'all','omitnan');
rMap = sum(spkmap,1,'omitnan');
pOcc = sum(bnoccs,1,'omitnan') / sum(bnoccs,'all','omitnan');
spatial_info = sum(pOcc .* rMap .* log2(rMap ./ uFR),1,'omitnan') ./ uFR;
end