function [root] = loadKS(datpath, spath, overwrite)
%% Loads the main file outputs from Phy
% Requires the npy-matlab package from kwikteam https://github.com/kwikteam/npy-matlab
% Generates or loads a root file in the ksPath and saves it to the starting
% directory
%
% Inputs:
% datpath   = string specifying the path to the kilosort output directory
% spath     = string specifying path where _root file will be saved
% overwrite = 0 or 1 to overwrite an existing '*_root' file in ksPath dir
%
% Outputs:
% root = struct containing information about the spike times and labels
%   fs      = sampling rate (30kHz default)
%   ts      = Mx1 array of all spike times in seconds
%   cl      = Mx1 array of cluster IDs per spike
%   lb      = Nx2 matrix of cluster IDs and labels from rater 
%   good    = Index of clusters assigned 'good' by rater
%   mua     = Index of clusters assigned 'mua' (Multi-Unit Acitivty) by rater
%   noise   = Index of clusters assigned 'noise' by rater
%
% Created 5/10/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    datpath
    spath
    overwrite = 0
end

cd(spath)

% Check if root file already exists and whether to overwrite
tmp = dir('*_root.mat');
if ~isempty(tmp) & overwrite == 0
    root = load(tmp.name); root = root.root;
    disp(['Loaded existing root file ', tmp.name])
    return
end

%% Get sync pulse by reading in meta data
cd(datpath)
lfpath = dir("*rec*");
cd(lfpath.name)

try
    syncfile = dir('*lf.bin');
    meta = SGLX_readMeta.ReadMeta(syncfile(1).name, syncfile(1).folder);
catch
    syncfile = dir('*tcat.imec0.ap.bin');
    meta = SGLX_readMeta.ReadMeta(syncfile(1).name, syncfile(1).folder); 
end

nSamp = floor(str2double(meta.fileTimeSecs) * SGLX_readMeta.SampRate(meta));
dataArray = SGLX_readMeta.ReadBin(0, nSamp, meta, syncfile(1).name, syncfile(1).folder);

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
% For 3B2 imec data: the sync pulse is stored in line 6.
dLineList = 6;

syncpulse = SGLX_readMeta.ExtractDigital(dataArray, meta, dw, dLineList);
syncpulse = double(syncpulse);
tspulse   = (1:length(syncpulse)) / SGLX_readMeta.SampRate(meta);

%% Get LFP data
if length(syncfile) > 1 %i.e. for 2.0 recordings
    chSubsamp = 4;
    lfSubsamp = 2;
    meta      = SGLX_readMeta.ReadMeta(syncfile(2).name, syncfile(2).folder);
    chanCt    = SGLX_readMeta.ChannelCountsIM(meta);    %[AP LFP SY] Chan counts
    dataArray = SGLX_readMeta.ReadBin(0, nSamp, meta, syncfile(2).name, syncfile(2).folder);
    dataArray = dataArray(1:chSubsamp:chanCt,1:lfSubsamp:length(dataArray)); %Subsample to 1250KHz and 1/chSubsamp channels
    dataArray = SGLX_readMeta.GainCorrectIM(dataArray, 1:size(dataArray,1), meta).*1000; %Outputs channels in mVolts
else
    chSubsamp = 4;
    lfSubsamp = 2;
    [~,chanCt]= SGLX_readMeta.ChannelCountsIM(meta);    %[AP LFP SY] Chan counts
    dataArray = dataArray(1:chSubsamp:chanCt,1:lfSubsamp:length(dataArray)); %Subsample to 1250kHz and 1/chSubsamp channels
    dataArray = SGLX_readMeta.GainCorrectIM(dataArray, 1:size(dataArray,1), meta).*1000; %Outputs channels in mVolts    
end

if contains(meta.imDatFx_pn, 'NP2_FLEX_0')
    meta.prbType = 'NPX1.0';
    root.prbType = 'NPX1.0';
elseif contains(meta.imDatFx_pn, 'NPM_FLEX_01')
    meta.prbType = 'NPX2.0';
    root.prbType = 'NPX2.0';
end

%% Get kilosort and Phy outputs
cd(datpath)
kspath = dir('*kilosort*');
try 
    cd([kspath.name '\sorter_output'])
catch
    cd(kspath.name)
end

fileClust = dir("spike_clusters.npy");
fileTime  = dir("spike_times.npy");
fileLabel = dir("cluster_group.tsv");   % Contains post-manual curation labels
fileInfo  = dir("cluster_info.tsv");    % Contains pre-manual curation labels and metadata like depth

try
    spkClusts   = readNPY(fileClust.name);
    spkTimes    = readNPY(fileTime.name);
    spkLabels   = readtable(fileLabel.name, "FileType", "text", 'Delimiter', '\t');
    spkInfo     = readtable(fileInfo.name, "FileType", "text", 'Delimiter', '\t');
catch
    disp("Missing critical file 'spike_clusters.npy', 'spike_times.npy', 'cluster_group.tsv', or 'cluster_info.tsv'. Aborting.")
    return
end

%% Organize root struct
[~, mname] = fileparts(meta.fileName);
mname = strsplit(mname, '_g0');
root.name = mname{1};
root.fs         = 30000;
root.ts         = double(spkTimes)/root.fs;
root.cl         = spkClusts;
root.lb         = spkLabels;
root.info       = spkInfo(:,[1:3,5:8,10]);
root.good       = root.lb.cluster_id(find(strcmp(root.lb.group,'good')));
root.mua        = root.lb.cluster_id(find(strcmp(root.lb.group,'mua')));
root.noise      = root.lb.cluster_id(find(strcmp(root.lb.group,'noise')));
root.goodind    = strcmp(root.lb.group,'good');
root.muaind     = strcmp(root.lb.group,'mua');
root.noiseind   = strcmp(root.lb.group,'noise');
root.syncpulse  = syncpulse;
root.tspulse    = tspulse;
root.fspulse    = SGLX_readMeta.SampRate(meta);

%% Assign shankID to units
if contains(meta.prbType,'NPX2.0')
    for i = 1:height(root.info)
        if root.info.ch(i) >= 0 && root.info.ch(i) < 48 || root.info.ch(i) >= 96 && root.info.ch(i) < 144
            root.info.shankID(i) = 0;
        elseif root.info.ch(i) >= 48 && root.info.ch(i) < 96 || root.info.ch(i) >= 144 && root.info.ch(i) < 192
            root.info.shankID(i) = 1;
        elseif root.info.ch(i) >= 192 && root.info.ch(i) < 240 || root.info.ch(i) >= 288 && root.info.ch(i) < 336
            root.info.shankID(i) = 2;
        elseif root.info.ch(i) >= 240 && root.info.ch(i) < 288 || root.info.ch(i) >= 336 && root.info.ch(i) < 384
            root.info.shankID(i) = 3;
        end
    end
    % root.info.shankID = floor(root.info.ch/96);
    % root.info.depth2 = root.info.depth - 720 .* root.info.shankID;
else 
    root.info.shankID = zeros(length(root.info.ch),1);
end

%% Process LFP and assign LFP metadata
% Use 3rd order Butterworth 300Hz lowpass zero phase lag filter
root.fs_lfp = SGLX_readMeta.SampRate(meta)/lfSubsamp;
[filtb, filta] = butter(3, 300/(root.fs_lfp/2),'low');
root.lfp       = filtfilt(filtb, filta, dataArray(1:size(dataArray,1),:)')';

lfpch = (1:chSubsamp:chanCt)'-1;   % Correct for 0-starting IDs of electrodes
lfpShank = zeros(size(root.lfp,1),1);
lfpDepth = zeros(size(root.lfp,1),1);
root.lfpinfo = table(lfpch,lfpShank,lfpDepth);

if contains(meta.prbType,'NPX2.0')
    for i = 1:height(root.lfpinfo)
        if root.lfpinfo.lfpch(i) >= 0 && root.lfpinfo.lfpch(i) < 48 || root.lfpinfo.lfpch(i) >= 96 && root.lfpinfo.lfpch(i) < 144
            root.lfpinfo.lfpShank(i) = 0;
        elseif root.lfpinfo.lfpch(i) >= 48 && root.lfpinfo.lfpch(i) < 96 || root.lfpinfo.lfpch(i) >= 144 && root.lfpinfo.lfpch(i) < 192
            root.lfpinfo.lfpShank(i) = 1;
        elseif root.lfpinfo.lfpch(i) >= 192 && root.lfpinfo.lfpch(i) < 240 || root.lfpinfo.lfpch(i) >= 288 && root.lfpinfo.lfpch(i) < 336
            root.lfpinfo.lfpShank(i) = 2;
        elseif root.lfpinfo.lfpch(i) >= 240 && root.lfpinfo.lfpch(i) < 288 || root.lfpinfo.lfpch(i) >= 336 && root.lfpinfo.lfpch(i) < 384
            root.lfpinfo.lfpShank(i) = 3;
        end
    end
    tmpDepthMap = probeDepthMap(root);  %Create depth map with potential shifts due to recording bank hot swapping
    root.lfpinfo.lfpDepth = tmpDepthMap(root.lfpinfo.lfpch+1)';
else 
    root.lfpinfo.lfpShank = zeros(length(root.lfpinfo.lfpch),1);
    tmpDepthMap = probeDepthMap(root);  %Create depth map with potential shifts due to recording bank hot swapping
    root.lfpinfo.lfpDepth = tmpDepthMap(root.lfpinfo.lfpch+1)';
end

%% Save output
cd(spath)

save([root.name, '_root'],'root')
save([root.name, '_meta'],'meta')

end

