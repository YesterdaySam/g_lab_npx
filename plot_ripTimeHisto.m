function [fhandle] = plot_ripTimeHisto(root,sess,bnsz)
%% Plots the histogram of appearance over a recording per shank
%
% Inputs:
% root = root object. Must have root.ripStruc field of struct(s) per shank
% sess = sess object
% bnsz = in sec, size of histogram bins
%
% Outputs:
% fhandle1 = handle to main IRI figure
% fhandle2 = handle to zoomed in 5 seconds IRI figure
%
% Created 5/31/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    sess
    bnsz = 30
end

nshanks = length(root.ripStruc);
cmapcool = cool(nshanks);
binedges = 0:bnsz:sess.ts(end)+bnsz;

fhandle = figure;
hold on

% Ripple appearance
for i = 1:nshanks
    plot(0,0,'Color',cmapcool(i,:))
    histogram(sess.ts(root.ripStruc(i).ripples(:,2)),binedges,'EdgeColor',cmapcool(i,:),'LineWidth',2,'DisplayStyle','stairs','HandleVisibility','off')
    legCell(i) = {['Sh' num2str(i-1)]};
end
xlabel('Recording Time (s)')
ylabel('Counts')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

end