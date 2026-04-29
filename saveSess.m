function [] = saveSess(sess,fname,spath)
% Saves a root file, takes in 'root' and optionally 'spath'

arguments
    sess
    fname = sess.name(1:14)
    spath = pwd;
end

sname = [spath, '\', fname, '_session'];
save(sname,'sess','-v7.3')

end