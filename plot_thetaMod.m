function [circ_stats,fhandle] = plot_thetaMod(root,unit,lfpInd,thbnsz,plotFlag)
%% Calculates polar rate relative to theta (0 = trough, 180 = peak)
% Methods from Oliva et al., 2016 (https://www.sciencedirect.com/science/article/pii/S0896627316305001#sec4)
% Requires Toolbox for circular statistics with Matlab. Authors: Philipp Berens; Email: philipp@bethgelab.org; Homepage: http://philippberens.wordpress.com/code/circstats/
%
% Inputs:
% root = root object. Must have root.thEnv field (bandpassed, hilbert
%   transformed theta)
% unit = cluster_id
% lfpInd = index of shankID e.g. 1 or 4
% plotFlag = binary. Whether to plot.
%
% Outputs:
% circ_stats = struct with the following fields:
%   mrl = mean resultant length of polar firing rate
%   ang = preferred firing direction (radian), 0 = trough, pi = peak
%   p = p-value of Rayleigh's test of non-uniformity for circular data
%   z = z score of Rayleigh's test
% fhandle = handle to figure
%
% Created 5/2/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    unit
    lfpInd % (use center of s.pyr. layer)
    thbnsz = 2*pi/18 % Default 20 degree bins (360/18 = 20)
    plotFlag = 1
end

spktsb = root.tsb(root.cl == unit);

% Get phase of analytic signal in theta range
% Originally, Theta peaks = 0/360 and troughs = 180/540
% thetaLFP = bandpass(root.lfp(lfpInd,:),root.bands(1,:), root.fs_lfp);
% th_env = abs(hilbert(thetaLFP));

% Offset theta by pi, for ease of calculations 0 to 2*pi
th_phase = angle(root.thEnv(lfpInd,:))+pi;
phase_edges = 0:thbnsz:2*pi;
th_occ = histcounts(th_phase,phase_edges)/root.fs_lfp;

% Assign spikes to lfp time bins
counts = histcounts(spktsb,root.lfp_tsb);
inds = [];

% find index of all bins >i spikes and concatenate to existing indices
for i = 0:max(counts)
    inds = cat(2,inds,find(counts>i));
end
% sort inds in case of repeats
inds = sort(inds);
inds = inds(:);

spk_phase = th_phase(inds);

spk_ct = histcounts(spk_phase,phase_edges);
polar_rate = spk_ct ./ th_occ;
scale_f = max(polar_rate); %ceil(max(polar_rate)/pi);

circ_ang = circ_mean(phase_edges(1:end-1)',polar_rate');
circ_mrl = circ_r(phase_edges(1:end-1)',polar_rate');
circ_mag = circ_mrl*scale_f;

[circ_p,circ_z] = circ_rtest(phase_edges(1:end-1)',polar_rate');

if plotFlag
    fhandle = figure;
    polarplot(phase_edges, [polar_rate polar_rate(1)],'k')
    hold on
    polarplot([circ_ang circ_ang], [0 circ_mag], 'r', 'LineWidth',2);
    scale_f = 3*pi/4*floor(scale_f/(2*pi));
    text(scale_f,scale_f, ['Rayleighs p = ' num2str(round(circ_p,3))], 'FontSize', 8);
end

% Prepare for output
circ_stats.mrl = circ_mrl;
circ_stats.ang = circ_ang;
circ_stats.p = circ_p;
circ_stats.z = circ_z;

% CMBHome code - duplicates circ toolbox, but shows math!
% xs = polar_rate.*cos(phase_edges(1:end-1)); % Mean Resultant length is sin and cosine
% ys = polar_rate.*sin(phase_edges(1:end-1));
% ang_hd = atan2(mean(ys),mean(xs)); % mean direction
% mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(polar_rate); % mean resultant length
% mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*2*pi;

% Plot raw lfp, filter lfp, theta envelope and phase, and unit spikes
% figure;
% plot(sess.ts(root.lfp_tsb),root.lfp(lfpInd,:),'k')
% hold on
% plot(sess.ts(root.lfp_tsb),thetaLFP,'r')
% plot(sess.ts(root.lfp_tsb),th_env,'g')
% plot(sess.ts(root.lfp_tsb),th_phase,'b')
% plot(sess.ts(spktsb),ones(length(spktsb),1)+3,'k|')

end