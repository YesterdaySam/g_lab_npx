function [root] = addRootLFP(spath)
%% Adds back extraneous LFP channels from lfp to root
%
% Inputs:
%   spath = string to directory with _root.mat and _lfp.mat files
%
% Outputs:
%   root = root object with original LFP data
%
% Created 1/15/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

rootF = dir('*_root.mat'); 
lfpF  = dir('*_lfp.mat'); 
load([spath '\' rootF.name])
load([spath '\' lfpF.name])

root.lfp     = lfp.lfp;
root.lfpinfo = lfp.lfpinfo;
root.uPSDMax = lfp.uPSDMax;
root.uPSD    = lfp.uPSD;

end