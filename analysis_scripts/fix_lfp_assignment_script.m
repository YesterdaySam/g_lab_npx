%% Fix LFP Channel assignment and layer boundaries

spath = 'D:\Data\Kelton\analyses\KW106\KW106_06252026_rec_D2_RLat1';
cd(spath)
rootfile = dir('*_root.mat');
load(rootfile.name)
sessfile = dir('*_session.mat');
load(sessfile.name)

ripRef = 4;
shanks = [4];
depths = [540];
lyrShs = [1,2,3,4];
lyrBds = [30  30  60 360; ...
          360 420 510 510];
% lyrShs = [1];
% lyrBds = [60; ...
%           240];

%% Add LFP channesl back to root and assign preferred LFP channels

root = addRootLFP(spath);

for i = 1:length(shanks)
    root = do_lfp2shank(root,shanks(i),depths(i));
end

%% Assign new layer boundaries

for i = 1:length(lyrShs)
    root.lyrbounds(:,lyrShs(i)) = lyrBds(:,i);
end

%% Plot new LFP peaks and bounds

lfpPowerFig = plot_lfpXdepth(root);
saveas(lfpPowerFig,[root.name '_LFPPower.png'])

%% Re assign units to layers and plot

root = get_layerUnits(root);

uTypeDepthFig = plot_layerUnits(root,1,0,0);
saveas(uTypeDepthFig,[root.name '_uTypeXDepth.png'])

%% Assign ripple reference shank and recalculate ripples

root.ripRef = ripRef;

chans = root.uPSDMax(2,:);
for i = 1:length(chans)
    ripStruc(i).ripples = get_ripples(root,chans(i),sess,3,5,[15 250]);
end
root.ripStruc = ripStruc;

%% Remove extra LFP channels and save
root = rmRootLFP(root);

saveRoot(root,spath)

disp('Fixed LFP and layer assignment')