function [binedges,binfr,spkmap] = get_frXpos(root,unit,sess,dbnsz,dend,vFlag,smFactor)
%% Collects smoothed firing rate by position and trial of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% dbnsz = size of position bins, default 0.05m = 5cm
% dend = length of track (m)
% vFlag = whether or not to remove spikes not coinciding with sess.runInds
% smFactor = how much to smooth the binned spikes and occupancy data
%
% Outputs:
% binedges = spatial bin edges
% binfr = MxN of M trials and N spatial-bins firing rates
% fhandle = handle to figure
%
% Created 7/23/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %m
    dend = 1.85
    vFlag = 1
    smFactor = 5;
end

% Only use valid trials
sess.lapstt = sess.lapstt(sess.valTrials);
sess.lapend = sess.lapend(sess.valTrials);
sess.nlaps  = length(sess.lapstt);
% sess.valTrials = 1:sess.nlaps;

binedges = 0:dbnsz:dend;
spkinds = root.tsb(root.cl == unit);

if vFlag
    spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
end

spkmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpspks = sess.pos(spkinds(spkinds > sess.lapstt(i) & spkinds < sess.lapend(i)));
    spkct   = histcounts(tmpspks, binedges);
    lapInds = sess.ind(sess.lapstt(i):sess.lapend(i));
    tmpRun  = lapInds(sess.runInds(lapInds));   % Use only run periods for calculating occupancy
    bnocc   = histcounts(sess.pos(tmpRun),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    spkmap  = [spkmap; spkct];              % Save spike counts
end

spksmooth = smoothdata(spkmap,2,'gaussian',smFactor);
occsmooth = smoothdata(bnoccs,2,'gaussian',smFactor);

binfr = spksmooth ./ occsmooth;
end