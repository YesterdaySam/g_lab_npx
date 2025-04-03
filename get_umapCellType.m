function [root] = get_umapCellType(root,umap_template,kmeans_centroids,c_IDs_IN, plotflag)
%% Calculates waveform width and FWHM and automatically assigns
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% umap_template = path to file containing umap template
% kmeans_centroids = path to file containing kmeans centroids
% c_IDs_IN = indices of kmeans centroids belonging to IN class(es)
% plotflag = binary of whether to plot the output
%
% Outputs:
% root = modified root object with root.info.uType putative unit type assignment
%       1 = Putative pyramidal
%       0 = Putative IN
%
% Created 4/2/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    umap_template = 'D:\Data\Kelton\analyses\group_analyses\waveform_classification\umap_euclidean_3.mat'
    kmeans_centroids = 'D:\Data\Kelton\analyses\group_analyses\waveform_classification\kmeans_centroids_3.mat'
    c_IDs_IN = [1]
    plotflag = 1
end

load('D:\Data\Kelton\analyses\group_analyses\waveform_classification\waveform_dat.mat')

% Align all unit peaks as negative
alignWFs = [];
for i = 1:size(root.templateWF,1)
    % Extract waveform peak
    tmpwf = root.templateWF(i,:);
    tmpmin = min(tmpwf);
    tmpmax = max(tmpwf);
    if abs(tmpmin) < tmpmax
        tmpwf = -tmpwf;
    end
    alignWFs = [alignWFs; tmpwf];
end

%Include pool of training waveforms to make clustering more robust
% wfIDs = logical([zeros(size(goodwfs,1),1); ones(size(alignWFs,1),1)]);
% poolWFs = [goodwfs; alignWFs];

% Run UMAP using template
reduction_wfs = run_umap(alignWFs, 'template_file', umap_template, 'verbose', 'none');
% reduction_wfs = reduction_wfs(wfIDs,:);

% K-means cluster based on old centroids 
load(kmeans_centroids);
[~,tmpidx] = pdist2(tmpc,reduction_wfs,'euclidean','Smallest',1);
% tmpidx = tmpidx(wfIDs);

% Assign spike clusters to IN clusters
INWFs = zeros(height(root.info),1);

for i = 1:length(c_IDs_IN)
    tmpINs = tmpidx == c_IDs_IN(i);
    INWFs(tmpINs) = 1;
end

%% Plot against FWHM ID'ed interneurons

if plotflag
    root = get_estCellType(root,15,5,100,0);

    cmapcool = cool(length(c_IDs_IN));

    umapINsFig = figure; hold on
    plot(reduction_wfs(root.info.uType == 0,1),reduction_wfs(root.info.uType == 0,2),'r*')
    plot(reduction_wfs(:,1),reduction_wfs(:,2),'k.')
    for i = 1:length(c_IDs_IN)
        plot(reduction_wfs(tmpidx == c_IDs_IN(i),1),reduction_wfs(tmpidx == c_IDs_IN(i),2),'.','Color',cmapcool(i,:))
    end
    xlabel('UMAP 1'); ylabel('UMAP 2'); xticks([]); yticks([]);
end

%% For testing purposes
% round_reds = round(reduction_wfs,3);
% tmpUnit = find(round_reds(:,1) == -5.881);
% figure; plot(root.templateWF(tmpUnit,:)); title(['Unit ' num2str(tmpUnit)]);

%% Return updated root

root.info.uType = ~INWFs;

end
