% Summarize ephys script

%% Load root from scratch
spath = 'D:\Data\Kelton\analyses\KW032\KW032_03062025_rec_D3_RMed1';
datpath = 'D:\Data\Kelton\probe_data\KW032\KW032_03062025_rec_D3_RMed1_g0';

loadKS(datpath,spath);
root = alignBhvrTS(spath,spath,spath);

%% Load existing root file
% spath = 'D:\Data\Kelton\analyses\KW031\KW031_02182025_rec_D4_LLat1';
saveFlag = 1;

cd(spath)

rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)

%%
% Plot Amp X Depth x FR
adfr1 = plot_ampXdepthxFR(root);
adfr2 = plot_ampXdepthxFR(root,1,1,0);
adfr3 = plot_ampXdepthxFR(root,1,0,0);
adfr4 = plot_ampXdepthxFR(root,0,1,0);

saveas(adfr1,[root.name '_AmpDepthFR_all.png'])
saveas(adfr2,[root.name '_AmpDepthFR_noNoise.png'])
saveas(adfr3,[root.name '_AmpDepthFR_good.png'])
saveas(adfr4,[root.name '_AmpDepthFR_mua.png'])

% Plot FR X Depth
frd1 = plot_frXdepth(root);
frd2 = plot_frXdepth(root,1,1,0);
frd3 = plot_frXdepth(root,1,0,0);
frd4 = plot_frXdepth(root,0,1,0);

if saveFlag
    saveas(frd1,[root.name '_FRDepth_all.png'])
    saveas(frd2,[root.name '_FRDepth_noNoise.png'])
    saveas(frd3,[root.name '_FRDepth_good.png'])
    saveas(frd4,[root.name '_FRDepth_mua.png'])
end

%% Summarize counts by shank and depth
tmpCountFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.1 0.6])
tmpedges = min(root.info.depth):20:max(root.info.depth);
binCounts = histcounts(root.info.depth(root.goodind),tmpedges);
barh(tmpedges(1:end-1)+10,binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410])
ylim([0 max(root.info.depth)])
xlabel('Good Counts')
ylabel('Distance from tip (um)');
set(gca,'FontSize',12,'FontName','Arial')

saveas(tmpCountFig,[root.name '_count_good.png'])

tmpCountShankFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.15])
tmpedges = 0:1:length(unique(root.info.shankID));
binCounts = histcounts(root.info.shankID(root.goodind),tmpedges);
b = bar(tmpedges(1:end-1),binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410]);
set(get(b,'Parent'),'ydir','reverse')
ylim([0 max(binCounts)])
xlim([-0.5 max(root.info.shankID)+0.5])
xlabel('Shank ID')
ylabel('Good Counts');
set(gca,'FontSize',12,'FontName','Arial')
text(tmpedges(1:end-1),binCounts,num2str(binCounts'),'vert','bottom','horiz','center','FontSize',12); 

if saveFlag
    saveas(tmpCountShankFig,[root.name '_shankCount_good.png'])
end

%% Summarize spike counts across session

tmpSpikeDensityFig = plot_spikeDensity(root,60);

if saveFlag
    saveas(tmpSpikeDensityFig,[root.name '_spikeDensity.png'])
end

%% Raw LFP power by depth

NFFT = 10000;
[pxx,f] = pwelch(root.lfp',[],[],NFFT,root.fs_lfp);

bands = [6 10; 150 250];    % Also [0.01 300]
uPSD = zeros(size(bands,1),size(pxx,2));

for band = 1:size(bands,1)
    f_tmp = f > bands(band,1) & f < bands(band,2);  % Get frequencies in band
    uPSD(band,:) = mean(10*log10(pxx(f_tmp,:)),1);    %Get mean per electrode within band
end

% Locate electrode with max spectral power for each band, for each shank
root.bands = bands;
for sh = 0:3
    for band = 1:size(bands,1)
        [tmpMax, tmpInd] = max(uPSD(band,root.lfpinfo.lfpShank == sh));
        tmpD = root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh);
        tmpD = tmpD(tmpInd);
        root.uPSDMax(band,sh+1) = find(root.lfpinfo.lfpDepth == tmpD & root.lfpinfo.lfpShank == sh);
    end
end
root.uPSD = uPSD;
root = get_layerUnits(root,100);    % Assign units to layers, min layer width 100um, prominence cut off 0.25

%% Plot LFP by depth

cmap(1).map = cool(4);
cmap(2).map = hot(8);
cmap(3).map = gray(5);

tmpLFPDensityFig = figure; hold on
legcell = {};
nChXSh = height(root.lfpinfo)/numel(unique(root.lfpinfo.lfpShank));
grid on
for sh = 0:3
    for band = 1:size(bands,1)
        plot(uPSD(band,root.lfpinfo.lfpShank == sh),root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh),'Color',cmap(band).map(sh+1,:))
        legcell(band) = {[num2str(bands(band,1)) '-' num2str(bands(band,2)) 'Hz']};
    end
    for band = 1:size(bands,1)
        [tmpMax, tmpInd] = max(uPSD(band,root.lfpinfo.lfpShank == sh));
        tmpD = root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh);
        plot(tmpMax+1,tmpD(tmpInd),'<','Color',cmap(band).map(sh+1,:))
        if band == 2; plot([tmpMax tmpMax]+1,root.lyrbounds(:,sh+1),':_','Color',cmap(band).map(sh+1,:)); end
    end
end
xlabel('PSD (dB/Hz)')
ylabel('Distance from Probe tip (um)')
legend(legcell)
set(gca,'FontSize',12,'FontName','Arial')

%% Save updated root

if saveFlag
    saveas(tmpLFPDensityFig,[root.name '_LFPPower.png'])
    saveRoot(root,spath)
end

%% Plot phys compact

close all

plotPhysCompact(root,sess,spath,1)

disp(['Finished phys summary for ' root.name])














