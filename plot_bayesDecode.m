function [fhandle] = plot_bayesDecode(decodeI)

% Visualize 
fhandle = figure; hold on;
set(gcf,'Units','normalized','Position',[0 0.65 1 0.15])
plot(decodeI.newT,decodeI.rPos,'b')
plot(decodeI.newT,decodeI.dPos,'r.')
plot(decodeI.newT,decodeI.dErr,'g')
plot([0 max(decodeI.newT)],[0 0],'k--')
legend('Real','Decode','Error')
ylabel('Position (cm)'); xlabel('time (sec)')
set(gca,'FontSize',16,'FontName','Arial')

end
