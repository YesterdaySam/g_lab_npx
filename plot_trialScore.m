function [fhandle] = plot_trialScore(sess,nBlock)
%% Plot Trial Success score moving average smoothed by nBlocks
% Inputs
%   sess    = struct from importBhvr.m
%   nBlock  = N span of trials to average over
%
% Outputs
%   fhandle  = handle to figure
%
% Created 4/1/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sess
    nBlock = 10     % nTrials per block
end

lapscore = get_trialSuccess(sess);
smoothScore = smooth(lapscore, nBlock);

fhandle = figure; hold on
plot(sess.valTrials, lapscore,'k-*','LineWidth',1)
plot(sess.valTrials, smoothScore,'r-','LineWidth',2)
legend({'Trial',[num2str(nBlock) ' trial avg.']},'Location','best')
xlabel('Trial #'); ylabel('Operant Success')
ylim([-0.05 1.05])
set(gca,'FontSize',12,'FontName','Arial')

end