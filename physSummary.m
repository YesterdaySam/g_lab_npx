% Summarize ephys script

%% Load root from scratch
spath = 'D:\Data\Kelton\analyses\ZM006\ZM006_10142025_rec_D2_RLat2';
datpath = 'D:\Data\Kelton\probe_data\ZM006\ZM006_10142025_rec_D2_RLat2_g0';

loadKS(datpath,spath,1); 
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
% adfr2 = plot_ampXdepthxFR(root,1,1,0);
% adfr3 = plot_ampXdepthxFR(root,1,0,0);
% adfr4 = plot_ampXdepthxFR(root,0,1,0);

% Plot FR X Depth
% frd1 = plot_datXdepth(root,0);
% frd2 = plot_datXdepth(root,0,1,1,0);
frd3 = plot_datXdepth(root,0,1,0,0);
% frd4 = plot_datXdepth(root,0,0,1,0);

if saveFlag
    saveas(adfr1,[root.name '_AmpDepthFR_all.png'])
    % saveas(adfr2,[root.name '_AmpDepthFR_noNoise.png'])
    % saveas(adfr3,[root.name '_AmpDepthFR_good.png'])
    % saveas(adfr4,[root.name '_AmpDepthFR_mua.png'])
    % saveas(frd1,[root.name '_FRDepth_all.png'])
    % saveas(frd2,[root.name '_FRDepth_noNoise.png'])
    saveas(frd3,[root.name '_FRDepth_good.png'])
    % saveas(frd4,[root.name '_FRDepth_mua.png'])
end

%% Summarize spike counts across session

tmpSpikeDensityFig = plot_spikeDensity(root,60);

if saveFlag
    saveas(tmpSpikeDensityFig,[root.name '_spikeDensity.png'])
end

%% Estimate LFP power by depth

root = get_lfpXdepth(root,[6 10; 150 250]);

% !====! Manually set theEnv channels if desired, otherwise comment out !====!
% chIDs = [14];
% nshanks = numel(unique(root.lfpinfo.lfpShank));
% for i = 1:nshanks
%     thetaLFPSh = bandpass(root.lfp(chIDs(i),:),root.bands(1,:), root.fs_lfp);
%     root.thEnv(i,:) = hilbert(thetaLFPSh);  % Operates along rows
% end

root = get_layerBounds(root,100,0.25);    % Assign units to layers, min layer width 100um, prominence cut off 0.25

%% Plot PSD (not very informative)

% figure; hold on
% nshank = size(root.uPSDMax,2);
% ind1hz = find(f >= 1,1);
% ind300hz = find(f >= 300,1);
% for sh = 1:nshank
%     cmapcool = cool(nshank);
%     plot(f(ind1hz:ind300hz),pxx(ind1hz:ind300hz,root.uPSDMax(2,sh)),'Color',cmapcool(sh,:))
% end

%% Plot LFP by depth

tmpLFPDensityFig = plot_lfpXdepth(root);

if saveFlag; saveas(tmpLFPDensityFig,[root.name '_LFPPower.png']); end

%% Assign putative unit type (0 = IN; 1 = Principle) and FR Variance by time

[root, INsFig] = get_estCellType(root,0.5,0.4,100,1);    % FW = 15; FWHM = 5; FR = 100; Plotflag = 1
% [root, INsFig] = get_umapCellType(root,1);
root = get_FRVar(root,sess,0);

if saveFlag; saveas(INsFig,[root.name '_FW_class.png']); end

%% Assign and plot units by layer and type

root = get_layerUnits(root);

uTypeDepthFig = plot_layerUnits(root,1,0,0);
if saveFlag; saveas(uTypeDepthFig,[root.name '_uTypeXDepth.png']); end

%% Summarize counts by shank and depth

[tmpCountFig, tmpCountShankFig] = plot_count_shank(root);

if saveFlag
    saveas(tmpCountFig,[root.name '_count_good.png'])
    saveas(tmpCountShankFig,[root.name '_shankCount_good.png'])
end

%% Get ripples at s. pyr. layer as determined by ripple band peak

nshanks = numel(unique(root.lfpinfo.lfpShank));

chans = root.uPSDMax(2,:);

ct = 1;
catRips = [];
for chan = chans
    ripStruc(ct).ripples = get_ripples(root,chan,sess,3,5,[15 250]);
    % catRips = [catRips; ripStruc(ct).ripples(:,2), -1+ct+zeros(size(ripStruc(ct).ripples,1),1)];
    ct = ct+1;
end

root.ripStruc = ripStruc;

%% Save updated root and split out LFP for subsequent root save/load speed

if saveFlag
    lfp.name     = root.name;
    lfp.lfp      = root.lfp;
    lfp.fs_lfp   = root.fs_lfp;
    lfp.lfpinfo  = root.lfpinfo;
    lfp.lfp_tsb  = root.lfp_tsb;
    lfp.bands    = root.bands;
    lfp.uPSDMax  = root.uPSDMax;
    lfp.uPSD     = root.uPSD;

    nshanks = numel(unique(root.lfpinfo.lfpShank));

    root.lfp     = root.lfp(root.uPSDMax(2,:),:);
    root.uPSD    = root.uPSD(:,root.uPSDMax(2,:));
    root.lfpinfo = root.lfpinfo(root.uPSDMax(2,:),:);
    root.uPSDMax = [repmat(1:nshanks,size(root.bands,2),1)]; % Set uPSD to match updated LFP info

    saveRoot(root,spath)
    save([lfp.name '_lfp'],'lfp','-v7.3')
end

%% Plot phys compact

close all

plotPhysCompact(root,sess,spath,1)

disp(['Finished phys summary for ' root.name])
