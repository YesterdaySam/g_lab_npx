function [] = saveRoot(root,spath)
% Saves a root file, takes in 'root' and optionally 'spath'

arguments
    root
    spath = pwd;
end

sname = [spath, '\', root.name, '_root'];
save(sname,'root')

end