% Import behavior script

bhvrdir = 'C:\Users\cornu\Documents\Research\Data\KW004\D3';

sess = importBhvr(bhvrdir);
sbase = sess.name(1:14);
save([sbase, '_session'], 'sess')

fig_trialvel = plot_trialvel(sess,1);  % 1cm binsize
saveas(fig_trialvel,[sbase, '_trialvelocity'],'png')

[fig_lickraster, fig_lickavg] = plot_lickraster_rwd(sess, 5, 0.1);
saveas(fig_lickraster,[sbase, '_lickraster_rwdalign'],'png')
saveas(fig_lickavg,[sbase, '_lickaverage'],'png')

if bhvrdir(end-1:end) ~= "D1"
    fig_lickpos = plot_lickpos(sess);   % Don't run this for D1 - random acclimation (# laps ~= # rewards)
    saveas(fig_lickpos,[sbase, '_lickraster'],'png')
end

disp(['Finished basic import and visualization of ', sbase])