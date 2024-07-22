function [binedges,binfr,fhandle] = plot_frXpos(root,unit,sess,dbnsz,plotflag)
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
    dbnsz = 0.03    %m
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
nlaps = size(sess.lapstt,1);

% binedges = 0:dbnsz:max(sess.pos);
binedges = 0:dbnsz:1.9; %Hardcoded for known track length <190

spkvel = sess.velshft(root.tsb(root.cl == unit));
runspk = spkvel > 0.02; %Only use spikes greater than 2cm/s

spkpos = sess.pos(root.tsb(root.cl == unit));

dspk = histcounts(spkpos(runspk),binedges);
docc = histcounts(sess.pos,binedges)/sess.samprate;
binfr = dspk./docc;

sem = std(binfr,'omitnan')/sqrt(nlaps);
ciup = rmmissing(mean(binfr,1,'omitnan') + sem*1.96);
cidn = rmmissing(mean(binfr,1,'omitnan') - sem*1.96);

if plotflag
    fhandle = figure; hold on
    plot([binedges(1:end-1)]*100,binfr, 'k')
    patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')

    xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')

    if max(binfr,[],'all') < 10
        ylim([0 10])
    elseif max(binfr,[],'all') < 20
        ylim([0 20])
    end

    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
end

end