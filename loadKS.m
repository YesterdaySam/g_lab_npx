function [root] = loadKS(ksPath, overwrite, sname)
%% Loads the main file outputs from Phy
% Requires the npy-matlab package from kwikteam https://github.com/kwikteam/npy-matlab
% Generates or loads a root file in the ksPath and saves it to the starting
% directory
%
% Inputs:
% ksPath    = string specifying the path to the kilosort output directory
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

parentDir = pwd;

fs = 30000; % hardcoded for NPX default sample rate 30KHz

cd(ksPath)

% Check if root file already exists and whether to overwrite
tmp = dir('*_root.mat');
if ~isempty(tmp) & overwrite == 0
    root = load(tmp.name); root = root.root;
    disp(['Loaded existing root file ', tmp.name])
    return
end

fileClust = dir("spike_clusters.npy");
fileTime = dir("spike_times.npy");
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
root.fs = fs;
root.ts = double(spkTimes)/fs;
root.cl = spkClusts;
root.lb = spkLabels;
root.good = find(strcmp(root.lb.group,'good'));
root.mua  = find(strcmp(root.lb.group,'mua'));
root.noise = find(strcmp(root.lb.group,'noise'));

cd(parentDir)

save([root.name, '_root'],'root')
end