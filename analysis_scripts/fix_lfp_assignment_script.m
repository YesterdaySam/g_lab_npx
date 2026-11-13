%% Fix LFP Channel assignment and layer boundaries

clear
close all

spath = 'D:\Data\Kelton\analyses\KW122\KW122_11052026_rec_D2_RMed1';

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
ripRef = 4; % 1-indexed
region = 1; % 1 = Sub/CA1 2 = EC
ec5stt = [];
ec3stt = []; % Use empty unless specifically editing
shanks = [4]; % Only use these for hand editing ref LFP depth per shank
depths = [570]; % Depth for hand editing of shank X
lyrShs = [1,2,3,4];
lyrBds = [0    60   210   360; ...
          270  480  510   600];

cd(spath)
rootfile = dir('*_root.mat');
load(rootfile.name)
sessfile = dir('*_session.mat');
load(sessfile.name)


%% Add LFP channels back to root and assign preferred LFP channels

root = addRootLFP(spath);

for i = 1:length(shanks)
    root = do_lfp2shank(root,shanks(i),depths(i));
end

%% Assign new layer boundaries

if region == 1
    for i = 1:length(lyrShs)
        root.lyrbounds(:,lyrShs(i)) = lyrBds(:,i);
    end
elseif region == 2
    [root,tmpDat] = get_lyrBounds_ec(root,sess,ec5stt,ec3stt);
end

%% Plot new LFP peaks and bounds

lfpPowerFig = plot_lfpXdepth(root);
fsave(lfpPowerFig,[root.name '_LFPPower'],1,0)

%% Re assign units to layers and plot

if region == 1
    root = get_layerUnits(root);
end

uTypeDepthFig = plot_layerUnits(root);
fsave(uTypeDepthFig,[root.name '_uTypeXDepth'],1,0)

if region == 2
    tmpD = root.info.depth(root.goodind);

    thDepthFig = plot_datXdepth(root,tmpDat,1,0,0,2); % theta
    plot(tmpDat.thAng(tmpDat.lyrID == 5)+180, tmpD(tmpDat.lyrID == 5),'k*')
    plot(tmpDat.thAng(tmpDat.lyrID == 3)+180, tmpD(tmpDat.lyrID == 3),'c*')
    plot(tmpDat.thAng(tmpDat.lyrID == 2)+180, tmpD(tmpDat.lyrID == 2),'m*')
    legend('off'); % legend('Data','EC5','EC3','EC2')
    set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]);
    thDepthFig = fixRatio(thDepthFig);
    fsave(thDepthFig,[root.name '_thXdepthXlyr'],1,0)
end

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