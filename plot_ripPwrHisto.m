function [fhandle] = plot_ripPwrHisto(root,bnsz,maxpwr)
%% Plots the histogram of ripple peak powers per shank
%
% Inputs:
% root = root object. Must have root.ripStruc field of struct(s) per shank
% bnsz = in db, size of histogram bins
% wlen = in db, maximum histogram bin
%
% Outputs:
% fhandle = handle to figure
%
% Created 5/31/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    bnsz = 0.025
    maxpwr = 0.75
end

nshanks = length(root.ripStruc);
cmapcool = cool(nshanks);

fhandle = figure;
hold on

binedges = 0:bnsz:maxpwr;

% Ripple peak power
for i = 1:nshanks
    plot(0,0,'Color',cmapcool(i,:))
    histogram(root.ripStruc(i).ripples(:,4),binedges,'EdgeColor',cmapcool(i,:),'LineWidth',2,'DisplayStyle','stairs','HandleVisibility','off')
    legCell(i) = {['Sh' num2str(i-1)]};
end

xlabel('Ripple Peak Power Amplitude')
xlim([0 maxpwr])
ylabel('Counts')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

end