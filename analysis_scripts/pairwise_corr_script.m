
% To ensure root.thEnv exists
root = get_lfpXdepth(root,[6 10; 150 250]);

%%
mainCC = 361;
mainSh = find(root.info.cluster_id == mainCC);
mainSh = root.info.shankID(mainSh) +1;

% plot_thetaMod(root,mainUnit,mainSh)

compIDs = find(root.goodind & root.info.lyrID == 1);

for i = 1:length(compIDs)
    close all

    compCC = root.info.cluster_id(compIDs(i));

    comp_acg = plot_acg(root,compCC);
    title(['Unit ' num2str(compCC)])
    comp_acg.Children.Children.FaceColor = [1 0 0];
    comp_ccg = plot_ccg(root,mainCC,compCC);
    comp_ccg.Children.Children(2).FaceColor = [0.5 0 0.5];
    main_acg = plot_acg(root,mainCC);
    title(['Unit ' num2str(mainCC)])

    figlist = get(groot, 'Children');
    newfig = figure;
    set(gcf,'units','normalized','position',[0.1 0.1 0.7 0.7])

    tcl = tiledlayout(newfig, 'flow');

    for j = 1:numel(figlist)
        figure(figlist(j))
        ax = gca;
        ax.Parent = tcl;
        ax.Layout.Tile = j;
    end
    figure(newfig)

    pause(2.5)
end