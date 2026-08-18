function [si,uFR,peakFR,peakLoc,spksmooth,occsmooth,binfr,binedges] = get_SI(root,unit,sess,dbnsz,dend)
%% Returns the Spatial Information of a Unit
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   unit = cluster ID
%   sess = session struct from importBhvr
%   dbnsz = size of position bins, default 0.05m = 5cm
%   dend = double in length meters of track length
%
% Outputs:
%   si = spatial info
%   uFR = mean FR over all bins
%   peakFR = peak FR
%   peakLoc = bin of FR peak
%   spksmooth = vector of circularly smoothed spike counts
%   occsmooth = vector of circularly smoothed occupancy time
%   binfr = vector of circularly smoothed FR
%   binedges = edges of spatial bins
%
% Created 10/31/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %m
    dend = 1.85     %m
end

binedges = 0:dbnsz:dend; %max(sess.pos(sess.lapstt(2):sess.lapend(2)));    % Base max binsize on first valid trial
nbins = length(binedges)-1;

spkinds = root.tsb(root.cl == unit);

try 
    spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods
catch
    disp('uh oh, failed to threshold spikes based on runInds')
    return
end

valspks = spkinds(sess.lapInclude(spkinds));
valoccs = sess.lapInclude & sess.runInds;

bnspks  = histcounts(sess.pos(valspks), binedges);
bnoccs  = histcounts(sess.pos(valoccs),binedges) / sess.samprate;

% Circular smooth
bnspks = repmat(bnspks,[1,3]);
bnoccs = repmat(bnoccs,[1,3]);

spksmooth = smoothdata(bnspks,'gaussian',5);
occsmooth = smoothdata(bnoccs,'gaussian',5);

spksmooth = spksmooth(nbins+1:nbins*2);
occsmooth = occsmooth(nbins+1:nbins*2);

binfr = spksmooth ./ occsmooth;
[peakFR,peakLoc] = max(binfr);

pOcc = occsmooth ./ sum(occsmooth,'all','omitnan');
uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
si = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan');

end