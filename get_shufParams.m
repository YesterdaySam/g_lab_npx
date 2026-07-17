function [datStruc] = get_shufParams(root,ccs,sess,datStruc,nShufs,doSpace,doSPWR,doCueFR,qStruc,histoBnsz,wlen)
%% Gets shuffled parameters for each unit in ccs based on flags
%
% Inputs:
%   root    = root struct
%   ccs     = vector of unit indices, such as root.good
%   sess    = session struct
%   datStruc= data struct
%   nShufs  = # of shuffles
%   doSpace = Spatial info, peak FR, peak loc, avg spatial FR, rate map and
%       spike map per trial
%   doSPWR  = average firing rate relative to ripple peak times
%   doCueFR = Cue info, cue avg FR, cue-relative rate map per trial
%   qStruc  = q onset/offset struc
%   dbnsz   = bin size for spatial calculations
%
% Outputs:
% datStruc = updated datStruc with shuffle parameters
%
% Created 5/26/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    ccs
    sess
    datStruc
    nShufs  = 250
    doSpace = true
    doSPWR  = false
    doCueFR = false
    qStruc  = 0
    histoBnsz = 5;
    wlen = 150;
end

tic

if root.ripRef < 1
    disp('Missing ripRef for this recording, skipping SPWR processing')
    doSPWR = false;
end

if doSpace
    datStruc.shufSI = zeros(length(ccs),nShufs);
end
if doSPWR
    datStruc.shufSPWR = zeros(nShufs,length(ccs),length(-wlen:histoBnsz:wlen)-1);
end
if doCueFR
    datStruc.shufQSI = zeros(length(ccs),nShufs);
end

for j = 1:nShufs
    [shiftroot,shiftsess] = shiftTrain(root,sess);

    if mod(j,50) == 0
        disp(['Shuffle # ' num2str(j)])
        toc 
    end

    for i = 1:length(ccs)
        cc = ccs(i);

        if doSpace
            datStruc.shufSI(i,j) = get_SI(shiftroot,cc,shiftsess);
        end

        if doSPWR
            datStruc.shufSPWR(j,i,:) = plot_frXripple(shiftroot,cc,shiftsess,root.ripRef,wlen,histoBnsz,0);
        end

        if doCueFR
            [datStruc.shufQSI(i,j)] = get_SI_cue(shiftroot,cc,shiftsess,qStruc);
        end
    end
end

toc

if doSpace
    datStruc.sigSI = sum(datStruc.shufSI > datStruc.trueSI,2) / nShufs;
end
if doCueFR
    datStruc.sigQSI = sum(datStruc.shufQSI > datStruc.qSI,2) / nShufs;
end


end