function [fhandle] = plot_templateWF(root,cc)
%% Plots the template waveform of a unit
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% unit = cluster ID
%
% Outputs:
% fhandle = handle to figure
%
% Created 2/12/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

ccc = find(root.info.cluster_id == cc);

maxY = max(root.templateWF(ccc,:));
minY = min(root.templateWF(ccc,:));

fhandle = figure; hold on
ts = (0:1:size(root.templateWF,2)-1) / root.fs * 1000;
plot(ts, root.templateWF(ccc,:),'k')
plot([ts(1) ts(end)], [0 0], 'k--')
ylim([minY-2, maxY+2])
ylabel('Amplitude (A.U.)')
xlabel('Time (ms)')
if root.info.uType(ccc) == 1
    typeStr = 'Putative Pyr';
else
    typeStr = 'Putative IN';
end
legend({['Probe Channel: ' num2str(root.info.ch(ccc))],typeStr})
set(gca,'FontSize',12,'FontName','Arial')

end