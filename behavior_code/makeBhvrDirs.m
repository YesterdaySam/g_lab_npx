function [] = makeBhvrDirs(mainDir)
%% Creates subdirectories for each behavior day and move files appropriately
% Sends the .h5 wavesurfer files and the .png files from Matlab to subdirs
%
% Inputs:
% mainDir = full path to a directory containing all behavior files
%
% Created 7/11/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

parentDir = pwd;

cd(mainDir)

imlist = dir('*.png');
wslist = dir('*.h5');

for i = 1:length(wslist)
    tmpfile = wslist(i).name;
    splitstr = tmpfile(1:end-8);
    mkdir(splitstr)
    newpath = [mainDir '\' splitstr];
    movefile(tmpfile, newpath)
end

dirlist = dir;

for i = 1:length(imlist)
    tmpim = imlist(i).name;
    splitstr = split(tmpim,'_laps.png');
    newpath = [mainDir '\' splitstr{1}];
    movefile(tmpim, newpath)
end

cd(parentDir)
end