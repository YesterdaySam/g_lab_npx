function [si,peakFR,uFR,pfields,binfr,binedges] = get_PF(root,unit,sess,dbnsz,vthresh)
%% Returns the SI, peak rate, and Place Field(s) of a Unit
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
% Created 11/01/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    dbnsz = 0.05    %m
    vthresh = 0.04  %m/s; velocity threshold for spikes
end

pfsizethresh = 0.15;

% nValLaps = length(sess.valTrials);

binedges = 0:dbnsz:max(sess.pos(sess.lapstt(2):sess.lapend(2)));    % Base max binsize on first valid trial
spkinds = root.tsb(root.cl == unit);
% spkinds = spkinds(sess.velshft(spkinds) > vthresh);     % Use only spikes above velocity threshold
spkinds = spkinds(sess.runInds(spkinds));   % Use only spikes in run periods

valspks = spkinds(sess.lapInclude(spkinds));
valoccs = sess.lapInclude & sess.runInds;

bnspks  = histcounts(sess.pos(valspks), binedges);

if sum(bnspks) < 100
    disp(['Fewer than 100 spikes found for unit ' num2str(unit) ' during get_PF'])
    si      = NaN;
    peakFR  = NaN;
    uFR     = NaN;
    pfields = NaN;
    binfr   = NaN;
    return
end

% bnoccs  = histcounts(sess.pos(sess.lapInclude),binedges) / sess.samprate;
bnoccs  = histcounts(sess.pos(valoccs),binedges) / sess.samprate;

spksmooth = smoothdata(bnspks,'gaussian',5);
occsmooth = smoothdata(bnoccs,'gaussian',5);

binfr = spksmooth ./ occsmooth;
peakFR = max(binfr);

% Find candidate fields based on 50% of peak smoothed FR
% fieldBins = binfr > 0.5*peakFR;

% Find candidate fields based on 3 std of mean
% fieldBins = binfr > mean + 3*std(binfr);

% Find candidate fields based on 3 std of noise floor (lower 50th pctile)
noiseLevel = binfr < prctile(binfr,50);
fieldBins = binfr > mean(binfr(noiseLevel)) + 3*std(binfr(noiseLevel));

fOffsets = find(diff([fieldBins 0]) == -1); 
fOnsets  = find(diff([0 double(fieldBins)]) == 1);
nFields = numel(fOffsets);

% For each candidate field, calculate basic metrics
for i = 1:nFields
    pfields(i).field = zeros(1,length(binfr));
    pfields(i).field(fOnsets(i):fOffsets(i)) = 1;
    pfields(i).field = logical(pfields(i).field);
    pfields(i).size = sum(pfields(i).field) * dbnsz;
    pfields(i).peak = max(binfr(pfields(i).field));
    pfields(i).infieldFR = mean(binfr(pfields(i).field));
    pfields(i).outfieldFR = mean(binfr(~pfields(i).field));
end

%Condense fields at very start and end of lap into 1 field for circularity
if pfields(1).field(1) == 1 && pfields(end).field(end) == 1
    pfields(1).field = logical(pfields(1).field + pfields(end).field);
    pfields(1).size = sum(pfields(1).field) * dbnsz;
    pfields(1).peak = max(binfr(pfields(1).field));
    pfields(1).infieldFR = mean(binfr(pfields(1).field));
    pfields(1).outfieldFR = mean(binfr(~pfields(1).field));
    pfields(end) = [];  %remove duplicate field
end

% Find fields < thresh width and remove
rmind = [];
for i = 1:length(pfields)
    if pfields(i).size < pfsizethresh
        rmind = [rmind i];
    end
end
pfields(rmind) = [];

% Get SI
pOcc = bnoccs ./ sum(bnoccs,'all','omitnan');
uFR = sum(spksmooth,'all','omitnan') / sum(occsmooth,'all','omitnan');
si = sum(pOcc .* binfr .* log2(binfr ./ uFR),'all','omitnan');

end