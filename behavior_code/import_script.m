% Batched Import behavior script

% bhvrdir = 'C:\Users\cornu\Documents\Research\Data\test\test_02';
% mousedir = 'C:\Users\cornu\Documents\Research\Data\KW005';
mousedir = 'D:\Kelton\ws_data\KW005';

overwriteFlag = 0;

cd(mousedir)

makeBhvrDirs(mousedir)

dirlist = dir;

for i = length(dirlist)-2:length(dirlist)
    cd(fullfile(dirlist(i).folder,dirlist(i).name));

    sessmade = dir('*session.mat');
    if ~isempty(sessmade) & overwriteFlag == 0
        disp(['Not a valid directory or session already created for ' dirlist(i).name])
        cd(mousedir)
        continue
    end

    try
        sess = importBhvr(fullfile(dirlist(i).folder,dirlist(i).name));
        sess = getErrorTrials(sess);
    catch
        cd(mousedir)
        continue
    end

    sbase = sess.name(1:14);
    save([sbase, '_session'], 'sess','-v7.3')
    disp(['Session created and saved for ' sbase])

    if isempty(sess.valTrials)
        continue
    end

    [fig_trialvel,fig_velavg, tmptrackedges, tmpbnvel] = plot_trialvel(sess,1);  % 1cm binsize
    saveas(fig_trialvel,[sbase, '_trialvelocity'],'png')
    saveas(fig_velavg,[sbase, '_avgvelocity'],'png')

    [fig_lickraster, fig_lickavg] = plot_lickraster_rwd(sess, 5, 0.1);
    saveas(fig_lickraster,[sbase, '_lickraster_rwdalign'],'png')
    saveas(fig_lickavg,[sbase, '_lickaverage'],'png')

    if dirlist(i).name(end-1:end) ~= "D1"
        [fig_lickpos, fig_licktrialavg, tmptrackedges2, ~, tmpbnlck] = plot_lickpos(sess);   % Don't run this for D1 - random acclimation (# laps ~= # rewards)
        saveas(fig_lickpos,[sbase, '_lickraster'],'png')
        saveas(fig_licktrialavg,[sbase, '_lickpos_average'],'png')
    end

    fig_vel_lick = plot_vel_lck(sess,tmpbnvel,tmpbnlck,tmptrackedges,tmptrackedges2);
    saveas(fig_vel_lick,[sbase, '_vel_lck'],'png')
    
    close all
    disp(['Finished basic import and visualization of ', sbase])

end

disp(['Finished basic import and visualization of all files for ' mousedir(end-4:end)])