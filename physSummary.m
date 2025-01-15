% Summarize ephys script

datpath = 'D:\Data\Kelton\analyses\KW022\KW022_12202024_rec_D4_LMed1';

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


%%
tmpCountFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.1 0.6])
tmpedges = min(root.info.depth):20:max(root.info.depth);
binCounts = histcounts(root.info.depth(root.goodind),tmpedges);
barh(tmpedges(1:end-1)+10,binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410])
ylim([0 max(root.info.depth)])
xlabel('Good Counts')
ylabel('Distance from tip (um)');
set(gca,'FontSize',12,'FontName','Arial')

saveas(tmpCountFig,[root.name '_count_good.png'])

tmpCountShankFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.15])
tmpedges = 0:1:length(unique(root.info.shankID));
binCounts = histcounts(root.info.shankID(root.goodind),tmpedges);
b = bar(tmpedges(1:end-1),binCounts,'FaceColor',[0.5 0.7235 0.8705],'EdgeColor',[0 0.4470 0.7410]);
set(get(b,'Parent'),'ydir','reverse')
ylim([0 max(binCounts)])
xlim([-0.5 max(root.info.shankID)+0.5])
xlabel('Shank ID')
ylabel('Good Counts');
set(gca,'FontSize',12,'FontName','Arial')
text(tmpedges(1:end-1),binCounts,num2str(binCounts'),'vert','bottom','horiz','center','FontSize',12); 

saveas(tmpCountShankFig,[root.name '_shankCount_good.png'])

disp(['Finished phys summary for ' root.name])

