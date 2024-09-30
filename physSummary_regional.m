% Summarize ephys script with more manual control over depth 

datpath = 'D:\Data\Kelton\analyses\KW007\KW007_08262024_rec_D1_Sub';

cd(datpath)

rootfile = dir("*_root.mat");
load(rootfile.name)

%% Manual input

sflag = 1;

border1 = 200;
border2 = 1000;

[root.ccs_g_ctx,root.ccinds_g_ctx] = getIDbyDepth(root,border2+1,max(root.info.depth));
[root.ccs_g_sub,root.ccinds_g_sub] = getIDbyDepth(root,border1+1,border2);
[root.ccs_g_dg, root.ccinds_g_dg]  = getIDbyDepth(root,0,border1);

ndg  = length(root.ccinds_g_dg);
nsub = length(root.ccinds_g_sub);
nctx = length(root.ccinds_g_ctx);

cmapdg  = winter(ndg);
cmapsub = copper(nsub);
cmapctx = summer(nctx);

legcell = {'Dentate gyrus','Subiculum','Cortex'};

%% Run block and save
if contains(root.prbType, 'NPX1.0')
    % Plot Amp X Depth x FR
    adfr1 = plot_regional_ADFR(root,root.ccinds_g_dg,root.ccinds_g_sub,root.ccinds_g_ctx,cmapdg,cmapsub,cmapctx);

    xs = xlim;
    plot([0,xs(2)],[border1 border1], 'k--')
    plot([0,xs(2)],[border2 border2], 'k--')

    legend(legcell,'Location','ne')

elseif contains(root.prbType, 'NPX2.0')
    

    disp(['Finished phys summary for ' root.name])
else
    disp(['Probe type not recognized for ' root.name])
end

if sflag
    saveas(adfr1,[root.name '_regional_ADFR_good.png'])
    disp(['Finished phys summary for ' root.name])
end
