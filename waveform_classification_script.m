%% IN Classification script
parentDir = "D:\Data\Kelton\analyses\group_analyses\Subiculum_FR"; 
subDatT = import_xldat(parentDir,"dat_include.xlsx");
cd(parentDir) 

cleanInds = subDatT.include == 0;
subDatT(cleanInds,:) = [];  %Clean excluded sessions

saveFlag = 1;

%% Load in data

goodwfs = [];
sessIDs = [];
goodfrs = [];

for i = 1:height(subDatT)
    cd(subDatT.fpath{i})

    rootfile = dir("*_root.mat");
    load(rootfile.name)
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    disp(root.name)

    nShank = numel(unique(root.info.shankID));

    goodwfs = [goodwfs; root.templateWF(root.goodind,:)];
    goodfrs = [goodfrs; root.info.fr(root.goodind,:)];
    sessIDs = [sessIDs; i*ones(sum(root.goodind),1)];
end

for i = 1:size(goodwfs,1)
    tmpwf = goodwfs(i,:);
    tmpmin = min(tmpwf);
    tmpmax = max(tmpwf);
    if abs(tmpmin) < tmpmax
        tmpwf = -tmpwf;
    end
    goodwfs(i,:) = tmpwf;
end

 %%
% goodwfs = root.templateWF(root.goodind,:);

unitFWHM = [];
unitFW   = [];
for i = 1:size(goodwfs)

    % Extract waveform peak
    tmpwf = goodwfs(i,:);
    tmpmin = min(tmpwf);
    tmpmax = max(tmpwf);
    if abs(tmpmin) > tmpmax
        tmpwf = -tmpwf;
    end

    % findpeaks(tmpwf,'Annotate','extents');    % Optional plotting
    [tmppks,tmplocs,tmpwidth,tmpprom] = findpeaks(tmpwf);
    [~,maxind] = max(tmppks);
    % [~,tmpinds] = sort(tmpprom);
    firsthalf = tmpwf(1:tmplocs(maxind));
    lasthalf  = tmpwf(tmplocs(maxind):end);
    lowind(1) = find(firsthalf <= tmppks(maxind) - tmpprom(maxind),1,'last');
    lowind(2) = find(lasthalf <= tmppks(maxind) - tmpprom(maxind),1,'first');

    unitFWHM = [unitFWHM; tmpwidth(maxind)];
    unitFW   = [unitFW; (lowind(2) + length(firsthalf)) - lowind(1)];
end

fwThresh = 15; fwhmThresh = 5;
INWFs = unitFW <= fwThresh & unitFWHM <= fwhmThresh;
% INWFs = root.info.uType(root.goodind) == 0;
INWFs = INWFs_UMAP;

%% Plot aggregated Wfs
figure; set(gcf,"Position",[680 300 495 630])
subplot(3,1,1:2); hold on
plot(unitFW(INWFs),unitFWHM(INWFs),'ro')
plot(unitFW(~INWFs),unitFWHM(~INWFs),'bo')
% plot(unitFW(root.good(umap_INs)),unitFWHM(root.good(umap_INs)),'go')
plot([fwThresh, fwThresh],[min(unitFWHM) max(unitFWHM)],'k--')
plot([min(unitFW) max(unitFW)],[fwhmThresh, fwhmThresh],'k--')
xlim([0 max(unitFW)]); ylim([0 max(unitFWHM)])
ylabel('FWHM (A.U.)')
set(gca,'FontSize',12,'FontName','Arial')

subplot(3,1,3); hold on
binctrs1 = 0:1:max(unitFW);
bincount = histogram(unitFW,binctrs1,'DisplayStyle','stairs','LineWidth',2,'EdgeColor','b');
incount = histogram(unitFW(INWFs),binctrs1,'DisplayStyle','stairs','LineWidth',2,'EdgeColor','r');
% bar(binctrs1(2:end)-0.5,bincount)
xlim([0 max(unitFW)]); xlabel('FullWidth (A.U.)'); ylabel('Count')
set(gca,'FontSize',12,'FontName','Arial')

%% PCA

[coefs,score,latent,tsquared,explained] = pca(goodwfs,'centered',false);

figure; hold on;
plot(explained,'k-o')
plot(cumsum(explained),'k-o')

xlim([0 length(explained)])
xlabel('PC#')
ylabel('Explained variance')

figure; 
plot3(score(:,1),score(:,2),score(:,3),'k.')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')

hold on
plot3(score(INWFs,1),score(INWFs,2),score(INWFs,3),'r.')
title(replace(root.name,'_',' '))

figure; hold on
plot(score(:,1),score(:,2),'k.')
plot(score(INWFs,1),score(INWFs,2),'r.')
xlabel('PC1')
ylabel('PC2')

%%
UMAP_dims = 2;  % Graphs created for 2 or 3.

% neighbor: Controls local and global structure
% low neighbors, more local. high, more global
UMAP_neighbor = 30; % Blake: 30  Wenbo: 50  Gardner et al. 2022: 5  UMAP default: 15

% min_dist: Controls how tightly UMAP is allowed to pack points.
UMAP_minDist = 0.3; % Blake: 0.6  Wenbo: 0.6  Gardner et al. 2022: 0.05  UMAP default: 0.3

% metric: controls how distance is computed in the ambient space
UMAP_metric = 'cityblock'; % 'cosine' (Tang Shin Jadhav); 'cityblock' (Gardner et al 2022); 'euclidean' (default); 'correlation'

%% run UMAP
% [reduction_wfs,umap] = run_umap(goodwfs,'n_components', UMAP_dims, 'metric', UMAP_metric, 'n_neighbors', UMAP_neighbor, 'min_dist', UMAP_minDist, 'verbose', 'none');
[reduction_wfs,umap] = run_umap(goodwfs, 'template_file', 'umap_euclidean_3.mat', 'verbose', 'none');
% [reduction_wfs,umap,clusters_wfs] = run_umap(goodwfs,'template_file', 'umap_template_1.mat', 'cluster_output', 'graphic', 'cluster_detail', 'medium');

if UMAP_dims == 2
    figure; hold on
    plot(reduction_wfs(INWFs,1),reduction_wfs(INWFs,2),'r*')
    plot(reduction_wfs(:,1),reduction_wfs(:,2),'k.')
    [tmpidx,tmpc] = kmeans(reduction_wfs,7);
    plot(reduction_wfs(tmpidx == 1,1),reduction_wfs(tmpidx == 1,2),'b.')
    plot(reduction_wfs(tmpidx == 2,1),reduction_wfs(tmpidx == 2,2),'g.')
    plot(reduction_wfs(tmpidx == 3,1),reduction_wfs(tmpidx == 3,2),'r.')
    plot(reduction_wfs(tmpidx == 4,1),reduction_wfs(tmpidx == 4,2),'c.')
    plot(reduction_wfs(tmpidx == 5,1),reduction_wfs(tmpidx == 5,2),'m.')
    xlabel('UMAP 1'); ylabel('UMAP 2'); xticks([]); yticks([]);
    legend('FWHM/FW INs')
    % figure; hold on
    % % cappedFRs = goodfrs;
    % % tmp = goodfrs > 40;
    % % cappedFRs(tmp) = 40;
    % scatter(reduction_wfs(:,1),reduction_wfs(:,2),[],goodfrs,'.')
    % colormap hsv
    % colorbar
    xlabel('UMAP 1'); ylabel('UMAP 2')
else
    figure;
    plot3(reduction_wfs(:,1),reduction_wfs(:,2),reduction_wfs(:,3),'k.')
    hold on
    plot3(reduction_wfs(INWFs,1),reduction_wfs(INWFs,2),reduction_wfs(INWFs,3),'r.')

end

round_reds = round(reduction_wfs,3);

% INWFs_UMAP = tmpidx == 2;

%%
tmpUnit = find(round_reds(:,1) == -4.230);
figure; plot(goodwfs(tmpUnit,:)); title(['Unit ' num2str(tmpUnit)]);

%% Run UMAP on new data using old template
% Must save template variable name as umap only otherwise generates errors!
fTemplate = 'D:\Data\Kelton\analyses\group_analyses\waveform_classification\umap_euclidean_3.mat';
[reduction_new] = run_umap(goodwfs,'template_file',fTemplate, 'verbose','none');

load('D:\Data\Kelton\analyses\group_analyses\waveform_classification\kmeans_centroids_3.mat');
tmpidx = kmeans(reduction_new, [], 'Distance','cityblock', 'Start', tmpc);

figure; hold on
plot(reduction_new(:,1),reduction_new(:,2),'k.')
plot(reduction_new(tmpidx == 1,1),reduction_new(tmpidx == 1,2),'b.')
plot(reduction_new(tmpidx == 2,1),reduction_new(tmpidx == 2,2),'g.')
plot(reduction_new(tmpidx == 3,1),reduction_new(tmpidx == 3,2),'r.')
plot(reduction_new(tmpidx == 4,1),reduction_new(tmpidx == 4,2),'c.')
plot(reduction_new(tmpidx == 5,1),reduction_new(tmpidx == 5,2),'m.')
xlabel('UMAP 1'); ylabel('UMAP 2'); xticks([]); yticks([]);

round_reds = round(reduction_wfs,3);

%% Plot different waveforms overlayed
nINs = sum(tmpidx == 3);
ct = 0;
goodINWfs = goodwfs(INWFs,:);
goodPyWfs = goodwfs(~INWFs,:);
% goodINWfs = goodwfs(tmpidx == 3,:);
% goodInvWfs = goodwfs(tmpidx == 2,:);
% goodPyrWfs = goodwfs(tmpidx == 4,:);
for i = 1:nINs/16
    figure
    for j = 1:16
        subplot(4,4,j)
        hold on
        plot(goodINWfs(j+ct,:),'k');
        plot(goodPyWfs(j+ct,:),'b')
        % plot(goodInvWfs(j+ct,:),'b')
        % plot(goodPyrWfs(j+ct,:),'r')
        plot([1,size(goodwfs,2)],[0 0],'k--')
        if j == 1
            legend({"IN Wf","Pyr Wf"})
        end
    end
    ct = ct +16;
end

%%