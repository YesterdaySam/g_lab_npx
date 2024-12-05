function [si,uFR,peakFR,spksmooth,occsmooth,fhandle] = get_SI(root,unit,sess,dbnsz,vthresh)
%% Returns the Spatial Information of a Unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% dbnsz = size of position bins, default 0.05m = 5cm
% vthresh = threshold of behavioral velocity to throw out spikes, default 0.04 m/s
%
% Outputs:
% binedges = spatial bin edges
% binfr = spatial-binned firing rate
% fhandle = handle to figure
%
% Created 10/31/24 LKW; Grienberger Lab; Brandeis University
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
% spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold
spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods

spkmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    spkct   = histcounts(tmpspks, binedges);
    lapInds = sess.ind(sess.lapstt(i):sess.lapend(i));
    tmpRun  = lapInds(sess.runInds(lapInds));   % Use only run periods for calculating occupancy
    % bnocc   = histcounts(sess.pos(sess.ind(sess.lapstt(i):sess.lapend(i))),binedges) / sess.samprate;
    bnocc   = histcounts(sess.pos(tmpRun),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);

binfr = spksmooth ./ occsmooth;
peakFR = max(binfr);

pOcc = occct ./ sum(occct,'all','omitnan');
uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
si = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan');

end