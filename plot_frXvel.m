function [binedges,binfr,mdlparams,fhandle] = plot_frXvel(spks,unit,sess,vbnsz,plotflag)
%% Plots the binned firing rate by velocity of a unit
%
% Inputs:
% spks = time series in seconds of spikes. Must be synced to sess.ts
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
% Created 7/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    spks            %spike time series
    unit {double}   %Cluster ID
    sess            %session struct
    vbnsz = 0.02    %m/s
    plotflag = 1    %binary
end

binedges = 0:vbnsz:max(sess.velshft);

% for i = 1:length(spks)
%     % spkvel(i) = sess.velshft(find(sess.ts > spks(i),1));
%     [~,indtmp] = min(abs(sess.ts - spks(i)));   %Find nearest velocity ts
%     spkvel(i) = sess.velshft(indtmp);
% end

spkvel = sess.velshft(spks);

vspk = histcounts(spkvel,binedges);
vocc = histcounts(sess.velshft,binedges)/sess.samprate;
binfr = vspk./vocc;

mdl = fitlm(binedges(1:end-1), binfr);
ys = predict(mdl,binedges(1:end-1)');
[r,p] = corrcoef(binfr,ys);
r = r(2,1);
p = p(2,1);
b = mdl.Coefficients{2,1} / 100; 
mdlparams.r = r;
mdlparams.p = p;
mdlparams.b = b;
mdlparams.yint = ys(1);

if plotflag
    fhandle = figure; hold on
    plot([binedges(1:end-1)]*100,binfr, 'ko','MarkerFaceColor','k')
    plot([binedges(1:end-1)]*100,ys,'r','LineWidth',2)
    xlabel('Velocity (cm/s)'); ylabel('Firing Rate')
    title(['Unit ' num2str(unit)])
    
    ylims = ylim;
    xlims = xlim;

    text(xlims(2) - .9*diff(xlims), ylims(2)-.1*diff(ylims), ['R = ' num2str(r, 3)], 'FontSize', 12)
    text(xlims(2) - .9*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(p, 3)], 'FontSize', 12)
    text(xlims(2) - .9*diff(xlims), ylims(2)-.2*diff(ylims), ['slope = ' num2str(b, 3)], 'FontSize', 12)
    text(xlims(2) - .9*diff(xlims), ylims(2)-.25*diff(ylims), ['y-int = ' num2str(ys(1), 3)], 'FontSize', 12)
end