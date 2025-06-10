function [si,uFR,peakFR,peakLoc,spksmooth,occsmooth,binfr,binedges] = get_SI(root,unit,sess,dbnsz,vthresh)
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

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(2):sess.lapend(2)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
% spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold
try 
    spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
catch
    disp('uh oh, failed to threshold spikes based on runInds')
    return
end

valspks = spkinds(sess.lapInclude(spkinds));
valoccs = sess.lapInclude & sess.runInds;

bnspks  = histcounts(sess.pos(valspks), binedges);
% bnoccs  = histcounts(sess.pos(sess.lapInclude),binedges) / sess.samprate;
bnoccs  = histcounts(sess.pos(valoccs),binedges) / sess.samprate;

spksmooth = smoothdata(bnspks,'gaussian',5);
occsmooth = smoothdata(bnoccs,'gaussian',5);

binfr = spksmooth ./ occsmooth;
[peakFR,peakLoc] = max(binfr);

pOcc = bnoccs ./ sum(bnoccs,'all','omitnan');
uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
si = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan');

end