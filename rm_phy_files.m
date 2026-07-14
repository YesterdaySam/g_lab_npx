function [] = rm_phy_files(parentdir,ksfolder)
%% Iterates through mouse and recording directories in parentdir and removes heavy phy files
% Inputs:
%   parentdir = path to directory containing either folders of many
%       recorded mice or individual recordings for a single mouse, or a
%       single recording ending in '_gX' where X is the rec iteration
% Created 5/7/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    parentdir               % Directory of 1 mouse or many mice
    ksfolder = 'kilosort4'  % Kilosort output subfolder (e.g. kilosort25)
end

curdir = pwd;
cd(parentdir)

dirlist = dir(parentdir);

if parentdir(end-2:end-1) == '_g' % For individual recordings

    for i = 3:length(dirlist)   % skip curdir and parent dir
        try
            rmdir(fullfile(dirlist(i).folder, ksfolder, 'sorter_output/.phy'), 's')
        catch
        end
    end

else    % For mouse folders or directories containing many mouse folders

    for i = 3:length(dirlist)   % skip curdir and parent dir
        cd(fullfile(dirlist(i).folder,dirlist(i).name))
        disp(['Checking directory: ' dirlist(i).name])

        try % For specific mouse directories (i.e. only phy folders in KW001)
            rmdir(fullfile(dirlist(i).folder, dirlist(i).name, ksfolder, 'sorter_output/.phy'), 's')
            disp(['Removed  .phy from: ' dirlist(i).name])
        catch
            subdirlist = dir(fullfile(dirlist(i).folder, dirlist(i).name));
            for j = 3:length(subdirlist)
                disp(['Checking sub-dirct: ' subdirlist(j).name])
                try % More generally for all phy files across all mice in parentdir
                    rmdir(fullfile(subdirlist(j).folder, subdirlist(j).name, ksfolder, 'sorter_output/.phy'), 's')
                    disp(['Removed  .phy from: ' subdirlist(j).name])
                catch
                    continue
                end
            end
        end
    end
end

disp('Completed removing .phy folder and subfolders')
cd(curdir)

end