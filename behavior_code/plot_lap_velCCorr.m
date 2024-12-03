function [spVelMat,f1,f2] = plot_lap_velCCorr(sess,dbnsz,vbnsz,plotflag)
%% Plots the cross correlation of average velocity per spatial bin between each lap
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields% unit = cluster ID
% sess = session struct from importBhvr
% dbnsz = size of histogram bins, default 0.05m = 5cm
% vbnsz = size of histogram bins, default 0.02m/s = 2cm/s
% plotflag = binary of whether to plot the output
%
% Outputs:
% binedges = spatial bin edges
% binfr = spatial-binned firing rate
% fhandle = handle to figure
%
% Created 12/2/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess            %session struct
    dbnsz = 0.05    %spatial bin size in m
    vbnsz = 0.002   %velocity bin in m/s
    plotflag = 1    %binary
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
nbins = length(binedges)-1;

spVelMat = zeros(sess.nlaps,nbins);

for i = 1:sess.nlaps
    tmppos = sess.pos(sess.lapstt(i):sess.lapend(i));
    tmpvel = sess.velshft(sess.lapstt(i):sess.lapend(i));
    for j = 1:nbins
        tmpin = tmppos > binedges(j) & tmppos < binedges(j+1);
        spVelMat(i,j) = mean(tmpvel(tmpin));
    end
end
[rhoPos,pPos] = corr(spVelMat,'Rows','complete');  %Test correlation between columns of matrix (spatially binned Velocity)
[rhoLap,pLap] = corr(spVelMat','Rows','complete');  %Test correlation between columns of matrix (laps)

if plotflag
    f1 = figure; imagesc(rhoPos); axis square;
    xticks(1:dbnsz*100:nbins)
    xticklabels(0:25:180)
    yticks(1:dbnsz*100:nbins)
    yticklabels(0:25:180)
    xlabel('Linearized Position (cm)')
    cbar1 = colorbar; %clim([0 inf]);
    ylabel(cbar1,'Pearson Rho','FontSize',12,'Rotation',90)

    % figure;imagesc(pPos < 0.05)

    f2 = figure; imagesc(rhoLap); axis square;
    xlabel('Lap #')
    cbar2 = colorbar; %clim([0 inf]);
    ylabel(cbar2,'Pearson Rho','FontSize',12,'Rotation',90)

    % figure;imagesc(pLap < 0.05)
end
end