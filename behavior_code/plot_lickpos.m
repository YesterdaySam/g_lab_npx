function [fhandle] = plot_lickpos(session)
%% Create linearized lick posiition (binned by space)
% Inputs
% session = struct from importBhvr.m
% Outputs
% fhandle = handle to figure

lckmap  = [];

for i = 1:session.nlaps-1
    % Find index of licks for this lap
    tmplck = session.lckind(find(session.lckind > session.lapstt(i) & session.lckind < session.lapend(i+1)));
    
    % Create Nx2 of [lick times, trial]
    lckmap = [lckmap; session.pos(tmplck), i*ones(numel(tmplck),1)];
end

if length(session.rwdind) ~= length(session.lapstt)
    for i = 1:session.nlaps-1
        tmpind(i) = find(session.rwdind > session.lapstt(i), 1);
    end
else
    tmpind = 1:session.nlaps;
end

rwdmap = [session.pos(session.rwdind), tmpind'];

fhandle = figure;      % Positional Lick Raster
hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
plot(lckmap(:,1)*100,lckmap(:,2),'k|')
plot(rwdmap(:,1)*100,rwdmap(:,2),'r*')
xlabel('Position (cm)')
ylabel('Trial #')
set(gca,'FontSize',12,'FontName','Arial')
title([strtok(session.name,'_'), ' ', session.name(end-9:end-8)])

end