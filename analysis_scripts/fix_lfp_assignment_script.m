%% Fix LFP Channel assignment and layer boundaries

clear
close all

spath = 'D:\Data\Kelton\analyses\KW117\KW117_10162026_rec_D2_RMed1';

% For copying hand editing from old root files
% oldpath = 'D:\Data\Kelton\analyses\ZM035\ZM035_06112026_rec_D1_LLat1';
% cd(oldpath)
% rootfile = dir('*_root.mat');
% load(rootfile.name)

% ripRef = root.ripRef
% shanks = [ripRef];
% depths = [root.lfpinfo.lfpDepth(ripRef)];
% lyrShs = [1,2,3,4];
% lyrBds = root.lyrbounds;

% For fresh hand editing
ripRef = 4;
shanks = []; % Only use these for hand editing ref LFP depth per shank
depths = []; % Depth for hand editing of shank X
lyrShs = [1,2,3,4];
lyrBds = [0    90   120   300; ...
          300   300  330   450];

cd(spath)
rootfile = dir('*_root.mat');
load(rootfile.name)
sessfile = dir('*_session.mat');
load(sessfile.name)


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