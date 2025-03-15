% Garbage Code from rippleDetect_script mostly subsumed by physSummary
% updates

%% FFT/PSD estimate of lfp
% clnlfp = bandstop(root.lfp(chan,:),[58 62],root.fs_lfp);

sigL = size(root.lfp,2);
chan = find(root.lfpinfo.lfpch == 336);
chans = find(root.lfpinfo.lfpShank == 0)';
% chans = 1:height(root.lfpinfo);

% Y = fft(root.lfp(chan,:));
% P2 = abs(Y/sigL);
% P1 = P2(1:sigL/2+1);
% f = root.fs_lfp/sigL*(0:(sigL/2));
% plot(f,P1,"LineWidth",1) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)"); xlim([0 300])
% ylabel("|P1(f)|")

% NFFT = 2^nextpow2(sigL);
NFFT = 10000;
[pxx,f] = pwelch(root.lfp(chans,:)',[],[],NFFT,root.fs_lfp);
% [pxx2,f] = pwelch(root.lfp(2,:),[],[],NFFT,root.fs_lfp);

pwrbands = [6 10];
fIncl = f>pwrbands(1) & f<pwrbands(2);
for i = 1:size(pxx,2)
    mPwr(i) = mean(10*log10(pxx(fIncl,i)));
end

%%
cmapcool = cool(size(chans,2));

figure; hold on
for i = 1:length(chans)
    tmpch = chans(i);
    plot(f,smooth(10*log10(pxx(:,tmpch)),20),'Color',cmapcool(i,:),"LineWidth",1) 
end
% plot(f,smooth(10*log10(pxx1),200),"LineWidth",1) 
% plot(f,smooth(10*log10(pxx2),200),"LineWidth",1) 
xlabel('Frequency (Hz)'); xlim([0 300])
ylabel('PSD (dB/Hz)')

%% Plot LFP power in bands across space on all shanks
bands = [6 10; 150 250];
uPSD = zeros(size(bands,1),size(pxx,2));

for band = 1:size(bands,1)
    f_tmp = f > bands(band,1) & f < bands(band,2);  % Get frequencies in band
    uPSD(band,:) = mean(10*log10(pxx(f_tmp,:)),1);    %Get mean per electrode within band
end

cmapcool = cool(4);
cmaphot = hot(8);
figure; hold on
for sh = 0:3
    plot(uPSD(1,root.lfpinfo.lfpShank == sh),root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh),'Color',cmapcool(sh+1,:))
    plot(uPSD(2,root.lfpinfo.lfpShank == sh),root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh),'Color',cmaphot(sh+1,:))
end

%%
bestChan = 4;
nRips = size(ripStruc(bestChan).ripples,1);
nearRips = zeros(nRips,nChans);
for i = 1:nRips
    pkInd = ripStruc(bestChan).ripples(i,2);
    % nearRips(i,1) = pkInd;
    for j = 1:nChans
        [tmpMin, tmpMinInd] = min(abs(pkInd - ripStruc(j).ripples(:,2)));
        if tmpMin < 1*root.fs_lfp     % 1 second
            nearRips(i,j) = ripStruc(j).ripples(tmpMinInd,2);
        else
            nearRips(i,j) = NaN;
        end
    end
end

%%
bnsz = 2.5;
compChans = logical(1 - ([1:nChans] == bestChan));
rDist = nearRips(:,compChans) - nearRips(:,bestChan);

bins = -root.fs_lfp:bnsz:root.fs_lfp;
binnedRipCounts = histcounts(reshape(rDist,1,[]),bins);

figure; hold on
% histogram(rDist)
bar((bins(1:end-1)+0.5*bnsz)/root.fs_lfp,binnedRipCounts);

figure; hold on
% plot(rDist(:,1)/root.fs_lfp,1:nRips,'k|')
for i = 1:nChans
    if compChans(i)
        plot(rDist(:,i)/root.fs_lfp,1:nRips,'|')
    end
end
plot((bins(1:end-1)+0.5*bnsz)/root.fs_lfp,binnedRipCounts,'r')
xlim([-0.02 0.02])
xlabel(['Time to reference ripple chan ' num2str(chans(bestChan))])
ylabel('Ripple #')