function [fhandle1,fhandle2] = plot_ripIRIHisto(root,sess,bnsz)
%% Plots the histogram of inter-ripple intervals per shank
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
    bnsz = 5
end

nshanks = length(root.ripStruc);
cmapcool = cool(nshanks);

fhandle1 = figure;
hold on

% Inter-Ripple Interval figures
for i = 1:nshanks
    IRIs = diff(sess.ts(root.ripStruc(i).ripples(:,1)));
    maxIRI(i) = max(IRIs);
    plot(0,0,'Color',cmapcool(i,:))
    histogram(IRIs, 0:bnsz:maxIRI(i)+bnsz,'EdgeColor',cmapcool(i,:),'LineWidth',2,'DisplayStyle','stairs','HandleVisibility','off')
    legCell(i) = {['Sh' num2str(i-1)]};
end
xlabel('Inter-Ripple-Interval (s)')
ylabel('Counts')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

fhandle2 = figure;
set(gcf,'units','normalized','position',[0.45 0.5 0.2 0.3])
hold on
for i = 1:nshanks
    IRIs = diff(sess.ts(root.ripStruc(i).ripples(:,1)));
    maxIRI(i) = max(IRIs);
    plot(0,0,'Color',cmapcool(i,:))
    histogram(IRIs, 0:0.25:5,'EdgeColor',cmapcool(i,:),'LineWidth',2,'DisplayStyle','stairs','HandleVisibility','off')
    legCell(i) = {['Sh' num2str(i-1)]};
end
xlabel('Inter-Ripple-Interval (s)')
ylabel('Counts')
legend(legCell)
set(gca,'FontSize',12,'FontName','Arial')

end