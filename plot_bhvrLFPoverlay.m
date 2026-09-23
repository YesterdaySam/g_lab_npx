function [fhandle] = plot_bhvrLFPoverlay(root, sess, twin, useSh, doVel, doPos, doLFPRaw, doTheta, doSPWR)
% Plots overlaid traces specified by doX flags over the time window in twin
% Position normalized to sess.maxPos, velocity normalized to 99th percentile
%
% Inputs:
%   root = root struct with root.ripStruc
%   sess = sess struct
%   twin = 1x2 double of indices specifying the time range to plot
%   useSh = shank to draw LFP from, defaults to root.ripRef
%   doVel = binary, add velocity trace
%   doPos = binary, add position trace
%   doLFPRaw = binary, add raw LFP trace
%   doTheta = binary, add filtered theta trace
%   doSPWR = binary, add filtered SPWR trace
%
% Outputs:
%   fhandle = handle to figure
%
% Created 8/21/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    sess
    twin
    useSh = root.ripRef
    doVel = 1
    doPos = 0
    doLFPRaw = 1
    doTheta = 1
    doSPWR = 1
end

% if mod(diff(twin),2) == 1   % Handle odd sized windows
%     twin(2) = twin(2)+1;
% end

indsLFP = find(root.lfp_tsb >= twin(1),1):1:find(root.lfp_tsb >= twin(2),1);
indsBhv = twin(1):20:twin(2);
xsLFP = sess.ts(root.lfp_tsb(indsLFP)) - sess.ts(root.lfp_tsb(indsLFP(1)));
xsBhv = sess.ts(indsBhv) - sess.ts(indsBhv(1));

fhandle = figure; hold on; axis off
fhandle.Renderer = 'Painters';
set(gcf,'units','normalized','position',[0.3536 0.4231 0.25 0.204])
fhandle = fixRatio(fhandle);

if doVel
    velTrace = sess.velshft(indsBhv) ./ 45;
    plot(xsBhv, velTrace, 'r')
end

if doPos
    posTrace = sess.pos(indsBhv) ./ sess.maxPos;
    plot(xsBhv, posTrace, 'c')
end

if doLFPRaw
    rawTrace = root.lfp(useSh,indsLFP) ./ max(abs(root.lfp(useSh,indsLFP))) .*2;
    plot(xsLFP, rawTrace + 1, 'k')
end

if doTheta
    thtlf = bandpass(root.lfp(useSh,:), [6 10], root.fs_lfp);
    thtTrace = thtlf(indsLFP) ./ max(abs(thtlf(indsLFP))) ./2; % prctile(abs(thtlf),99.9); 
    plot(xsLFP, thtTrace + 2, 'g')
end

if doSPWR
    riplf = bandpass(root.lfp(useSh,:), [150, 250], root.fs_lfp);
    ripEnv = abs(hilbert(riplf));
    ripTrace = riplf(indsLFP) ./ max(abs(riplf(indsLFP))) ./2; %  prctile(abs(riplf),99.9);
    plot(xsLFP, ripTrace + 3, 'b')
    
    % ripEnvTrace = ripEnv(indsLFP) ./ prctile(abs(ripEnv),99.9); % max(abs(riplf(indsLFP))) ./2;
    % plot(xsLFP, ripEnvTrace + 3, 'r')
    
    % non-functional, not finding ripple peaks correctly
    % ripPks = find(root.ripStruc(useSh).ripples(:,2) > indsLFP(1) & root.ripStruc(useSh).ripples(:,2) < indsLFP(end));
    % ripPks = root.ripStruc(useSh).ripples(:,2);
    % plot(xsLFP(ripPks), ones(size(ripPks)), 'r*')
end

% Plot scale bar
plot([0 xsLFP(end)/4], [3.5 3.5], 'k')
plot([0 0], [2.5 3.5], 'k')

end