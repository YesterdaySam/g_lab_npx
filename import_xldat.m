function [datT] = import_xldat(fpath,fname)

parentDir = pwd;
cd(fpath)

datT = readtable(fname);

% datT = baseT(:,1:4);
% datT.fpath = baseT.fpath;
% datT.include = baseT.include;
% datT.mouse = baseT.mouse;
% datT.session = baseT.session;

cd(parentDir)
end