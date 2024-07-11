function [fhandle] = plot_lickpos2(session)
%% Create linearized lick posiition (binned by space)
% Warning - doesn't handle first-laps without rewards very well, offsets
% actual reward delivery from licks by X number of laps where X =
% unrewarded laps
% Inputs
% session = struct from importBhvr.m
% Outputs
% fhandle = handle to figure

lckmap  = [];

% Use only rewarded laps
for i = 1:session.nlaps-1
    tmpind = session.rwdind > session.lapstt(i) & session.rwdind < session.lapend(i+1);
    isrewarded(i) = logical(sum(tmpind));
end
session.lapstt = session.lapstt(isrewarded);
session.lapend = session.lapend(isrewarded);

% Find licks
for i = 1:session.nlaps-1
    % Find index of licks for this lap
    tmplck = session.lckind(find(session.lckind > session.lapstt(i) & session.lckind < session.lapend(i+1)));
    
    % Create Nx2 of [lick times, trial]
    lckmap = [lckmap; session.pos(tmplck), i*ones(numel(tmplck),1)];
end

for i = 1:session.
% Correct for any deliveries prior to first true lap start
if session.ts(session.rwdind(1)) < session.ts(session.lapstt(1))
    session.rwdind = session.rwdind(2:end);
end

% Correct for any unrewarded laps
for i = 1:session.nlaps-1
    tmpind = session.rwdind > session.lapstt(i) & session.rwdind < session.lapend(i+1);
    isrewarded(i) = logical(sum(tmpind));
end
session.lapstt = session.lapstt(isrewarded);
session.lapend = session.lapend(isrewarded);

% Correct for any remaining differences between rwd delivery and lap #s
if length(session.rwdind) ~= length(session.lapstt)
    for i = 1:session.nlaps-1
        tmpind(i) = find(session.rwdind > session.lapstt(i), 1);
    end
else
    tmpind = 1:length(session.lapstt);
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