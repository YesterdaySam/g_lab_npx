function [si,uFR,peakFR,spksmooth,occsmooth,fhandle] = get_SI(root,unit,sess,dbnsz,vthresh)
%% Plots the avg binned firing rate by position of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% dbnsz = size of position bins, default 0.05m = 5cm
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
    dbnsz = 0.05    %m
    vthresh = 0.04  %m/s; velocity threshold for spikes
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(1):sess.lapend(1)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold

spkmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    spkct   = histcounts(tmpspks, binedges);
    bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);

binfr = spksmooth ./ occsmooth;
peakFR = max(binfr);

% rawfr = spkct ./ occct;

% if plotflag
%     fhandle = figure; hold on
%     plot([binedges(1:end-1)]*100,rawfr,'r')
%     plot([binedges(1:end-1)]*100,mean(binfr,1,'omitnan'), 'k')
%     % patch(100*[binedges(1:length(cidn)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
% 
%     xlabel('Position (cm)'); ylabel('Firing Rate (spk/s)')
% 
%     if max(mean(binfr,1,'omitnan'),[],'all') < 10
%         ylim([0 10])
%     elseif max(binfr,[],'all') < 20
%         ylim([0 20])
%     end
% 
%     title(['Unit ' num2str(unit)])
%     set(gca,'FontSize',12,'FontName','Arial')
% 
% end

% uFR = sum(spkct,'all','omitnan') / sum(occct,'all','omitnan');
% rMap = sum(spkmap,1,'omitnan');
pOcc = occct ./ sum(occct,'all','omitnan');
% spatial_info = sum(pOcc .* rawfr .* log2(rawfr ./ uFR),'all','omitnan') ./ uFR;

uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
si = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan') ./ uFR;

