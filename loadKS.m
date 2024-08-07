function [root] = loadKS(datpath, spath, sname, overwrite)
%% Loads the main file outputs from Phy
% Requires the npy-matlab package from kwikteam https://github.com/kwikteam/npy-matlab
% Generates or loads a root file in the ksPath and saves it to the starting
% directory
%
% Inputs:
% ksPath    = string specifying the path to the kilosort output directory
% sname     = string specifying name e.g. 'KW004_06282024'
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
    sname
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

% Get sync pulse by reading in meta data
cd(datpath)

lfpfile = dir('*lf.bin');
meta = SGLX_readMeta.ReadMeta(lfpfile.name, lfpfile.folder);

nSamp = floor(str2double(meta.fileTimeSecs) * SGLX_readMeta.SampRate(meta));
dataArray = SGLX_readMeta.ReadBin(0, nSamp, meta, lfpfile.name, lfpfile.folder);

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
% For 3B2 imec data: the sync pulse is stored in line 6.
dLineList = 6;

syncpulse = SGLX_readMeta.ExtractDigital(dataArray, meta, dw, dLineList);
syncpulse = double(syncpulse);
tspulse    = [1:length(syncpulse)] / SGLX_readMeta.SampRate(meta);

% Get kilosort and Phy outputs
kspath = dir('*kilosort*');
cd(kspath.name)

fileClust = dir("spike_clusters.npy");
fileTime  = dir("spike_times.npy");
fileLabel = dir("cluster_group.tsv");

try
    spkClusts   = readNPY(fileClust.name);
    spkTimes    = readNPY(fileTime.name);
    spkLabels   = readtable(fileLabel.name, "FileType", "text", 'delimiter', '\t');
catch
    disp("Missing critical file 'spike_clusters.npy', 'spike_times.npy', or 'spike_times.npy'. Aborting.")
    return
end

% Organize root struct
root.name = sname;  % Change later to use auto-generated filenames from parent dir
root.fs         = 30000;
root.ts         = double(spkTimes)/root.fs;
root.cl         = spkClusts;
root.lb         = spkLabels;
root.good       = root.lb.cluster_id(find(strcmp(root.lb.group,'good')));
root.mua        = root.lb.cluster_id(find(strcmp(root.lb.group,'mua')));
root.noise      = root.lb.cluster_id(find(strcmp(root.lb.group,'noise')));
root.syncpulse  = syncpulse;
root.tspulse    = tspulse;
root.fspulse     = SGLX_readMeta.SampRate(meta);

cd(spath)

save([root.name, '_root'],'root')
end