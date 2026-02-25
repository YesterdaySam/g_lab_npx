function [binedges,binBR,bstmap,uBR] = get_burstsXpos(root,unit,sess,dbnsz,dend,vFlag,smFactor)
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
% binBR = MxN of M trials and N spatial-bins burst rates
% bstmap = MxN of M trials and N spatial bin burst counts
% uBR   = 1xN of mean burst rate over N spatial-bins
%
% Created 12/4/25 LKW; Grienberger Lab; Brandeis University
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
[bStts] = get_bursts(root,sess,unit);

if vFlag
    bStts = bStts(sess.runInds(bStts));   % Use only spikes in run periods
end

bstmap = [];
bnoccs = [];
for i = 1:sess.nlaps
    tmpbst = sess.pos(bStts(bStts > sess.lapstt(i) & bStts < sess.lapend(i)));
    bstct   = histcounts(tmpbst, binedges);
    lapInds = sess.ind(sess.lapstt(i):sess.lapend(i));
    tmpRun  = lapInds(sess.runInds(lapInds));   % Use only run periods for calculating occupancy
    bnocc   = histcounts(sess.pos(tmpRun),binedges) / sess.samprate;
    bnoccs  = [bnoccs; bnocc];              % Save bin occupancy
    bstmap  = [bstmap; bstct];              % Save burst counts
end

bstsmooth = smoothdata(bstmap,2,'gaussian',smFactor);
occsmooth = smoothdata(bnoccs,2,'gaussian',smFactor);

binBR = bstsmooth ./ occsmooth;
uBR = sum(bstsmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');

end