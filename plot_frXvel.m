function [binedges,binfr,mdlparams,fhandle] = plot_frXvel(root,unit,sess,vbnsz,plotflag)
%% Plots the binned firing rate by velocity of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% vbnsz = size of velocity bins, default 2cm/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% binedges = velocity bin edges
% binfr = velocity-binned firing rate
% mdlparams = R and p values of the correlation coefficient, slope and
%   y-intercept of a linear model fit between velocity and firing rate
% fhandle = handle to figure
%
% Created 7/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    vbnsz = 2       %cm/s
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);

binedges = 0:vbnsz:max(sess.velshft);

spkvel = sess.velshft(root.tsb(root.cl == unit));

vspk = histcounts(spkvel,binedges);
vocc = histcounts(sess.velshft,binedges)/sess.samprate;
binfr = vspk ./ vocc;
lowOcc = vocc < 1; % Ignore bins with <1sec occupancy
lowOcc(1:2) = logical([1 1]); % Exclude periods <4cm/s velocity
binfr(lowOcc) = [];

% usebins = binedges > prctile(binedges,5) & binedges < prctile (binedges, 95);
usebins = logical([0 1-lowOcc(1:end)]);  %Drop first 2 bins (4cm/s) of standing periods and low occ bins
binctrs = binedges(usebins) - vbnsz/2;
% binfr = binfr(usebins);

% mdl = fitlm(binedges(1:end-diff([length(binedges) length(binfr)])), binfr);
% ys = predict(mdl,binedges(1:end-diff([length(binedges) length(binfr)]))');
mdl = fitlm(binctrs, binfr);
ys = predict(mdl, binctrs');
[r,p] = corrcoef(binfr,ys);
r = r(2,1);
p = p(2,1);
b = mdl.Coefficients{2,1}; 
mdlparams.r = r;
mdlparams.p = p;
mdlparams.b = b;
mdlparams.yint = predict(mdl,0);

if plotflag
    fhandle = figure; hold on
    % plot([binedges(1:end-diff([length(binedges) length(binfr)]))],binfr, 'ko','MarkerFaceColor','k')
    % plot([binedges(1:end-diff([length(binedges) length(binfr)]))],ys,'r','LineWidth',2)
    plot(binctrs,binfr, 'ko','MarkerFaceColor','k')
    plot(binctrs,ys,'r','LineWidth',2)
    xlabel('Velocity (cm/s)'); ylabel('Firing Rate (spk/s)')
    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')

    ylims = ylim;
    ylim([0 ylims(2)])
    ylims = ylim;
    xlims = xlim;
    xlim([0 xlims(2)])
    xlims = xlim;

    if b > 0
        text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims), ['R = ' num2str(r, 3)], 'FontSize', 12)
        text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(p, 3)], 'FontSize', 12)
        text(xlims(2) - .9*diff(xlims), ylims(2)-.2*diff(ylims), ['slope = ' num2str(b, 3)], 'FontSize', 12)
        text(xlims(2) - .9*diff(xlims), ylims(2)-.25*diff(ylims), ['y-int = ' num2str(mdlparams.yint, 3)], 'FontSize', 12)
    else
        text(xlims(2) - .3*diff(xlims), ylims(2)-.1*diff(ylims), ['R = ' num2str(r, 3)], 'FontSize', 12)
        text(xlims(2) - .3*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(p, 3)], 'FontSize', 12)
        text(xlims(2) - .3*diff(xlims), ylims(2)-.2*diff(ylims), ['slope = ' num2str(b, 3)], 'FontSize', 12)
        text(xlims(2) - .3*diff(xlims), ylims(2)-.25*diff(ylims), ['y-int = ' num2str(mdlparams.yint, 3)], 'FontSize', 12)        
    end
end

end