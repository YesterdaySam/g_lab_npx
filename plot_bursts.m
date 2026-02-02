function [fhandle] = plot_bursts(root,sess,unit)
% Plot spike times of all spikes in all bursts for a unit
%
% Inputs
%   root = root struct
%   sess = session struct
%   unit = single unit ID
%
% Outputs:
%   fhandle   = handle to figure
%
% Created 12/4/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

[bStts,~,bLens] = get_bursts(root,sess,unit);

spkinds = root.tsb(root.cl == unit);
spkts = sess.ts(spkinds);

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.17 0.32])
cols = turbo(max(bLens));
for i = 1:length(bStts)
    burstAlign = find(spkinds == bStts(i));
    tmpBurst = spkts(burstAlign:burstAlign + bLens(i) -1);
    tmpAlign = tmpBurst - sess.ts(bStts(i));
    scatter(tmpAlign*1000,i+zeros(size(tmpAlign)),[],cols(1:bLens(i),:))
end
set(gca,'FontSize',12,'FontName','Arial')
xlabel('Time from burst onset (ms)')
xlim([-.001 inf])
ylabel('Burst #')

end