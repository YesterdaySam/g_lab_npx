function [] = fixvelshft_batch(mousedir)

cd(mousedir)

dirlist = dir;
dirlist = dirlist([dirlist.isdir]);

for i = 1:length(dirlist)
    cd(fullfile(dirlist(i).folder,dirlist(i).name));

    try
        sessfile = dir("*session.mat");
        load(sessfile.name)

        sess.velshft = sess.velshft * 80;   %converts old value to cm/s

        save([sess.name(1:14) '_session.mat'],'sess')
        disp(['Updated and saved velocity for session file in ' dirlist(i).name])
    catch
        disp(['Failed to update velshft for directory ' dirlist(i).name])
    end
end

end