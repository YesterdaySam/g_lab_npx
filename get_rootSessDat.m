function [root,sess,dat] = get_rootSessDat(sdir,getRoot,getSess,getDat)
%% Plot and save some neural analyses
%
% Inputs:
%   sdir    = directory with root, sess, or data files
%   getRoot = binary, whether to try to load the root
%   getSess = binary, whether to try to load the sess
%   getDat  = binary, whether to try to load the data structure
%
% Outputs:
%   root    = root data
%   sess    = sess data
%   dat     = data struct
%
% Created 10/14/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    sdir
    getRoot = 1
    getSess = 1
    getDat  = 1
end

cd(sdir)

if getRoot
    try
        rootfile = dir("*_root.mat");
        load(rootfile.name)
        disp(['Loaded ' rootfile.name])
    catch
        disp('Failed to load root')
    end
else
    root = [];
end

if getSess
    try
        sessfile = dir("*_session.mat");
        load(sessfile.name)
        disp(['Loaded ' sessfile.name])
    catch
        disp('Failed to load sess')
    end
else
    sess = [];
end

if getDat
    try
        datfile = dir("*_dat.mat");
        load(datfile.name)
        disp(['Loaded ' datfile.name])
    catch
        disp('Failed to load data struct')
    end
else
    dat = [];
end

end