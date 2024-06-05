function [root] = loadKS(ksPath)
%%% loads the main file outputs from Phy including spike_times.npy,
%%% spike_clusters.npy and cluster_group.tsv
% Created 5/10/24 LKW; Grienberger Lab; Brandeis University
% Requires the npy-matlab package from kwikteam https://github.com/kwikteam/npy-matlab
% Inputs:
% ksPath = string specifying the path to the kilosort output directory
% Outputs:
% root = struct containing information about the 

parentDir = pwd;

fs = 30000; % hardcoded for NPX default sample rate 30KHz

cd(ksPath)

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

root.ts = double(spkTimes)/fs;
root.cl = spkClusts;
root.lb = spkLabels;
root.good = strcmp(root.lb.group,'good');
root.mua  = strcmp(root.lb.group,'mua');
root.noise = strcmp(root.lb.group,'noise');

cd(parentDir)
end