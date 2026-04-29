function [cueSI,uFR,pkFR,pkLoc,spkmap,bnoccs,rastMt,binfr,binedges] = get_SI_cue(root,unit,sess,qstruc,dbnsz)
%% Returns the Spatial Information of a Unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
% sess = session struct from importBhvr
% qstruc = qstruc.qs, qstruc.stt, qstruc.end incides, qstruc.qID trial IDs
% dbnsz = size of position bins, default 0.05m = 5cm
%
% Outputs:
%
% Created 4/17/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    qstruc
    dbnsz = 0.05    %m
end

try
    beltLength = sess.maxPos;
catch
    beltLength = 1.85;
end
wlen = beltLength / 2; % meters
binedges = -wlen:dbnsz:wlen;    % 100 cm (80 on each side of rand cue)
spkinds = root.tsb(root.cl == unit);
spkinds = spkinds (sess.runInds(spkinds));   % Use only spikes in run periods

spkmap = [];
bnoccs = [];
rastMt = [];

for i = 1:length(qstruc.q)
    % Find spikes, align to rand cue and bin
    tmpspks = spkinds(spkinds > qstruc.stt(i) & spkinds < qstruc.end(i));
    pstLapSpks = tmpspks > sess.lapstt(qstruc.qID(i)+1);
    preLapSpks = tmpspks < sess.lapend(qstruc.qID(i)-1);
    tmpspkpos = sess.pos(tmpspks) - sess.pos(qstruc.q(i));
    tmpspkpos(pstLapSpks) = tmpspkpos(pstLapSpks) + beltLength;
    tmpspkpos(preLapSpks) = tmpspkpos(preLapSpks) - beltLength;

    % Calculate occupancy over run periods and align to rand cue
    lapInds = sess.ind(qstruc.stt(i):qstruc.end(i));
    tmpRun  = lapInds(sess.runInds(lapInds));
    pstLapRun = tmpRun > sess.lapstt(qstruc.qID(i)+1);
    preLapRun = tmpRun < sess.lapend(qstruc.qID(i)-1);
    tmpPos = sess.pos(tmpRun) - sess.pos(qstruc.q(i));
    tmpPos(pstLapRun) = tmpPos(pstLapRun) + beltLength;
    tmpPos(preLapRun) = tmpPos(preLapRun) - beltLength;

    % Store lap info
    bnoccs = [bnoccs; histcounts(tmpPos,binedges) / sess.samprate];
    spkmap = [spkmap; histcounts(tmpspkpos,binedges)];
    rastMt = [rastMt; tmpspkpos qstruc.qID(i)*ones(size(tmpspkpos))];
end

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

spksmooth = smoothdata(spkct,'gaussian',5);
occsmooth = smoothdata(occct,'gaussian',5);
ratemap   = spkmap ./ bnoccs;

%needs NaN to be 0
binfr = spksmooth ./ occsmooth;
binfr(isnan(binfr)) = 0;
[pkFR,pkLoc] = max(binfr);

% Calculate SI relative to cue
pOcc = occsmooth ./ sum(occsmooth,'all','omitnan');
uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
cueSI = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan');

end