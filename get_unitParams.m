function [datStruc] = get_unitParams(root,ccs,sess,datStruc,doSpace,doVel,doTheta,doSPWR,doCueFR,qStruc,doOpto,dbnsz,wlen,histoBnsz)
%% Gets parameters for each unit in ccs based on flags
%
% Inputs:
%   root    = root struct
%   ccs     = vector of unit indices, such as root.good
%   sess    = session struct
%   datStruct = existing data struct, otherwise creates a fresh one
%   doSpace = Spatial info, peak FR, peak loc, avg spatial FR, rate map and
%       spike map per trial
%   doVel   = velocity linear model
%   doTheta = spikes relative to theta statistics
%   doSPWR  = average firing rate relative to ripple peak times
%   doCueFR = Cue info, cue avg FR, cue-relative rate map per trial
%   qStruc  = q onset/offset struc
%   doOpto  = Get opto-relative FR per trial and first peak after pulses
%   dbnsz   = bin size for spatial calculations
%
% Outputs:
% datStruc = datStruc
%
% Created 5/21/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    ccs
    sess
    datStruc = struct;
    doSpace  = true
    doVel    = true
    doTheta  = true
    doSPWR   = true
    doCueFR  = false
    qStruc   = 0
    doOpto   = false
    dbnsz    = 0.05
    wlen     = 150
    histoBnsz= 5
end

datStruc.name   = root.name;

if root.ripRef < 1
    disp('Missing ripRef for this recording, skipping SPWR processing')
    doSPWR = false;
end

% Arrange vector sizes properly
if doSpace
    datStruc.trueSI = zeros(length(ccs),1);
    datStruc.truePk = zeros(length(ccs),1);
    datStruc.trueLc = zeros(length(ccs),1);
end
if doSPWR
    datStruc.swrfr  = zeros(length(ccs),length(-wlen:histoBnsz:wlen)-1);
end
if doCueFR
    datStruc.qSI    = zeros(length(ccs),1);
end

% Get Requested parameters
for i = 1:length(ccs)
    cc = ccs(i);

    if doSpace
        [datStruc.trueSI(i),~,datStruc.truePk(i),datStruc.trueLc(i),~,~,datStruc.posfr(i,:),datStruc.binedges] = get_SI(root,cc,sess,dbnsz);
        [~,datStruc.frMap(:,:,i),datStruc.spkMap(:,:,i)] = get_frXpos(root,cc,sess,0.05,sess.maxPos,1);
    end

    if doVel
        [~,~,datStruc.trueVelMdl(i)] = plot_frXvel(root,cc,sess,2,0);
    end

    if doTheta
        lfpInd = root.info.shankID(root.info.cluster_id == cc)+1;   % Account for 0-indexing
        [datStruc.thetastats(i),datStruc.thetafr(i,:)] = plot_thetaMod(root,cc,lfpInd,2*pi/36,0);
    end

    if doSPWR
        [datStruc.swrfr(i,:),~,datStruc.swrz(i,:)] = plot_frXripple(root,cc,sess,root.ripRef,wlen,histoBnsz,0);
    end

    if doCueFR
        [datStruc.qSI(i),~,~,~,datStruc.qFRMap(:,:,i),~,datStruc.qPosFR(i,:)] = get_SI_cue(root,cc,sess,qStruc);
    end

    if doOpto
        [datStruc.optobins,datStruc.optoMat(i,:,:)] = plot_frXopto(root,cc,sess,0.002,0.5,0);
        datStruc.optoPkT(i) = get_firstPk(datStruc,i);
    end
end

if doSPWR
    datStruc.ripRate = size(root.ripStruc(root.ripRef).ripples,1) / (sum(not(sess.runInds)) / sess.samprate);    %Normalize based on standing periods
end

datStruc.binpos = datStruc.binedges(1:end-1)+0.5*dbnsz;

end