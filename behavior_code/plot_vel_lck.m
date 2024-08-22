function [fhandle] = plot_vel_lck(sess,bnvel,bnszVel,bnszLck)
%% Create linearized lick posiition (binned by space)
% Inputs
%   sess    = struct from importBhvr.m
%   bnvel   = binned velocity per trial output from plot_trialvel.m
%   bnszVel = double in meters (m)
%   bnszLck = double in meters (m)
%
% Outputs
%   fhandle = handle to figure
%
% Created 8/21/2024 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

nLaps = size(bnvel,1);   %Recalculate based on laps used in bnvel, don't use sess.nlaps
tracklen    = (max(sess.pos) - min(sess.pos)); % m
edgesVel    = 0:bnszVel:tracklen;
edgesLck    = 0:bnszLck:tracklen;
% nBinVel     = length(edgesVel);
nBinlck     = length(edgesLck);

% Find licks
lckmap  = [];
bnlck = ones(nLaps,nBinlck-1);

for i = 1:nLaps
    % Find index of licks for this lap
    tmplck = sess.lckind(find(sess.lckind > sess.lapstt(i) & sess.lckind < sess.lapend(i)));
    
    % Create Nx2 of [lick times, trial]
    lckmap = [lckmap; sess.pos(tmplck), i*ones(numel(tmplck),1)];
    
    bnocc = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),edgesLck);
    bnlck(i,:) = histcounts(sess.pos(tmplck),edgesLck);
    bnlck(i,:) = bnlck(i,:) ./ (bnocc / sess.samprate);  %Normalize to time in each bin
end


semvel = std(bnvel,'omitnan')/sqrt(nLaps);
ciupvel = rmmissing(mean(bnvel,1,'omitnan') + semvel*1.96);
cidnvel = rmmissing(mean(bnvel,1,'omitnan') - semvel*1.96);

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
plot(edgesVel(1:end-1)*100,mean(bnvel,1,'omitnan'),'k','LineWidth',2)
patch(100*[edgesVel(1:length(cidnvel)),fliplr(edgesVel(1:length(cidnvel)))],[cidnvel,fliplr(ciupvel)],'k','FaceAlpha',0.5,'EdgeColor','none')
xlabel('Position'); xlim([0 200])
ylabel('Average Velocity'); ylim([0 prctile(sess.velshft,99)*100])
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

semlck = std(bnlck,'omitnan')/sqrt(nLaps);
ciuplck = rmmissing(mean(bnlck,1,'omitnan') + semlck*1.96);
cidnlck = rmmissing(mean(bnlck,1,'omitnan') - semlck*1.96);

yyaxis right
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.3])
plot(edgesLck(1:end-1)*100,mean(bnlck,1,'omitnan'),'b','LineWidth',2)
patch(100*[edgesLck(1:length(cidnlck)),fliplr(edgesLck(1:length(cidnlck)))],[cidnlck,fliplr(ciuplck)],'b','FaceAlpha',0.5,'EdgeColor','none')
ylabel('Average Licks/s'); ylim([0 inf])
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

tmp = gca;
tmp.YAxis(2).Color = 'b';
tmp.YAxis(1).Color = 'k';

end