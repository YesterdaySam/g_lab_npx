function [fhandle] = plot_ripDurHisto(root,sess,bnsz,wlen)
%% Plots the histogram of ripple durations per shank
%
% Inputs:
% root = root object. Must have root.ripStruc field of struct(s) per shank
% sess = sess object
% bnsz = in ms, size of histogram bins
% wlen = in ms, maximum histogram bin
%
% Outputs:
% fhandle = handle to figure
%
% Created 5/31/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    sess
    bnsz = 5
    wlen = 150
end

nshanks = length(root.ripStruc);
cmapcool = cool(nshanks);

fhandle = figure;
hold on

binedges = 0:bnsz:wlen;

% plot Ripple durations
for i = 1:nshanks
    ripStts = sess.ts(root.ripStruc(i).ripples(:,1));
    ripStps = sess.ts(root.ripStruc(i).ripples(:,3));
    ripDurs = ripStps - ripStts;
    plot(0,0,'Color',cmapcool(i,:))
    histogram(ripDurs*1000,binedges,'EdgeColor',cmapcool(i,:),'LineWidth',2,'DisplayStyle','stairs','HandleVisibility','off')
    legCell(i) = {['Sh' num2str(i-1)]};
end

xlim([0 wlen])
xlabel('Ripple Duration (ms)')
ylabel('Counts')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

end