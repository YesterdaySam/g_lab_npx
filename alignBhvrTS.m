function [root] = alignBhvrTS(bhvrdir, rootdir, sdir)
%% Aligns neural time stamps to behavior using flipper pulse
% Uses cross correlogram of sync pulses sampled by two different clocks
% (SGLX clock and Wavesurfer clock) to identify offset of the random signal
% Then aligns all edges and updates spike timestamps to the nearest flipper
% pulse edge in the behavior stream. Finally bins all newly aligned spikes
% into the nearest behavior timestamp and saves a new root file
%
% Inputs:
%   bhvrdir    = path to directory containing *_sess.mat file
%   rootdir    = path to directory containing *_root.mat file
%   sdir       = path to save root file
%
% Outputs:
%   root = updated root object with tssync and tsb variables added
%       tssync = Nx1 double of spike times aligned to behavior clock
%       tsb    = Nx1 double of spike times binned by behavior time series
%
% Created 8/21/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

cd(rootdir)
rootfile = dir('*_root*');
rootstruc = load(rootfile.name); root = rootstruc.root; clear rootstruc;

cd(bhvrdir)
sessfile = dir('*_sess*');
sessstruc = load(sessfile.name); sess = sessstruc.sess; clear sessstruc;

cd(sdir)

[sync_ccg, sync_lags] = xcorr(root.syncpulse, sess.slx);

pk_lag = sync_lags(sync_ccg == max(sync_ccg));  %Index to offset sess.ts

[~, edges_sync] = findpeaks(root.syncpulse);
edges_root_ts = [root.tspulse(1) root.tspulse(edges_sync)];   %Start at sync.ts wall time (not a proper down/up flip)
% tmp = diff(edges_root_ts);

[~, edges_sess] = findpeaks(sess.slx);  %Indices of edges
align_start = edges_sess + pk_lag > 0;
edges_sess = [abs(pk_lag) edges_sess(align_start)'];   %Indices of sess edges after aligning to start of SGLX rec
edges_sess_ts = sess.ts(edges_sess);

%% align all spike times to behavior
% SA = EA + RA*(SB - EB)/RB
% SA is scaled, bounded, timestamp of event from B mapped to stream A
% SB is timestamp of event from B
% E is TS of Nearest Rising Edge in A or B
% R is sampling rate of stream A or B

allts = root.ts;
root.tssync = allts;
ra = sess.samprate;
rb = root.fspulse;

for i = 1:length(allts)
    try
        eb_ind = find(edges_root_ts < allts(i),1,'last');
        eb = edges_root_ts(eb_ind);   % Root sync edge ts
        % ea = edges_sess_ts(find(edges_sess_ts > eb,1,'first'));         % First sess sync edge after eb
        ea = edges_sess_ts(eb_ind);
    catch
        eb = edges_root_ts(1);
        ea = edges_sess_ts(find(edges_sess_ts > eb,1,'first'));
        disp(['During spike alignment used first sync ts edge to align spike ' num2str(i)])
    end

    try
        root.tssync(i) = ea + ra*(allts(i) - eb)/rb;
    catch
        eb = edges_root_ts(1);
        ea = edges_sess_ts(find(edges_sess_ts > eb,1,'first'));
        root.tssync(i) = ea + ra*(allts(i) - eb)/rb;
        % root.tssync(i) = allts(i);
        disp(['During spike alignment used first sync edge ts for spike ' num2str(i)])
    end
end

%% Account for alignment causing spikes to fall outside sess.ts(end)
tsovershoot = find(root.tssync > sess.ts(end));
if ~isempty(tsovershoot)
    root.ts(tsovershoot) = [];
    root.cl(tsovershoot) = [];
    root.tssync(tsovershoot) = [];
    disp(['During spike alignment ' num2str(length(tsovershoot)) ' spikes fell after sess.ts(end) and were removed'])
    root.trunc_spks_n = length(tsovershoot);
else
    root.trunc_spks_n = 0;
end

%% split spikes into bins with edges of sess.ts

counts = histcounts(root.tssync,sess.ts);   %after aligning all spike ts use root.tssync
inds = [];

% find index of all bins >i spikes and concatenate to existing indices
for i = 0:max(counts)
    inds = cat(2,inds,find(counts>i));
end
% sort inds in case of repeats
inds = sort(inds);
inds = inds(:);

root.tsb = inds;    %tsb = behavioral time series spikes

% Naively align LFP to start of sess, using pk_lag only (unbounded error)
root.lfp_tsb = sess.ind+pk_lag;
root.lfp_tsb = sess.ind(root.lfp_tsb > 0);
root.lfp_tsb = root.lfp_tsb(1:length(root.tspulse)); 
if length(root.lfp_tsb) > length(sess.ind)
    disp(['During LFP alignment, ' num2str(length(root.lfp_tsb) - length(sess.ind)) ' indices fell after sess.ind(end) and were removed'])
    root.lfp_tsb = root.lfp_tsb(1:length(sess.ind));
end

save([root.name '_root'],'root');

end