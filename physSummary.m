% Summarize ephys script

datpath = 'D:\Data\Kelton\analyses\KW007\KW007_09042024_rec_D6_Lat2';

cd(datpath)

rootfile = dir("*_root.mat");
load(rootfile.name)

% Plot Amp X Depth x FR
adfr1 = plot_ampXdepthxFR(root);
adfr2 = plot_ampXdepthxFR(root,1,1,0);
adfr3 = plot_ampXdepthxFR(root,1,0,0);
adfr4 = plot_ampXdepthxFR(root,0,1,0);

saveas(adfr1,[root.name '_AmpDepthFR_all.png'])
saveas(adfr2,[root.name '_AmpDepthFR_noNoise.png'])
saveas(adfr3,[root.name '_AmpDepthFR_good.png'])
saveas(adfr4,[root.name '_AmpDepthFR_mua.png'])

% Plot FR X Depth
frd1 = plot_frXdepth(root);
frd2 = plot_frXdepth(root,1,1,0);
frd3 = plot_frXdepth(root,1,0,0);
frd4 = plot_frXdepth(root,0,1,0);

saveas(frd1,[root.name '_FRDepth_all.png'])
saveas(frd2,[root.name '_FRDepth_noNoise.png'])
saveas(frd3,[root.name '_FRDepth_good.png'])
saveas(frd4,[root.name '_FRDepth_mua.png'])

disp(['Finished phys summary for ' root.name])
