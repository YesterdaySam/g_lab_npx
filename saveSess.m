function [] = saveSess(sess,spath)
% Saves a root file, takes in 'root' and optionally 'spath'

arguments
    sess
    spath = pwd;
end

sname = [spath, '\', sess.name(1:14), '_session'];
save(sname,'sess','-v7.3')

end