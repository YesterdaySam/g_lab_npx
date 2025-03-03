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
plot(root.templateWF(ccc,:),'k')
plot([1 length(root.templateWF(ccc,:))], [0 0], 'k--')
ylim([minY-2, maxY+2])
ylabel('Amplitude (A.U.)')
legend(['Probe Channel: ' num2str(root.info.ch(ccc))])
set(gca,'FontSize',12,'FontName','Arial')

end