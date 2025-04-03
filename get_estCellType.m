function [root] = get_estCellType(root,fwThresh,fwhmThresh,frThresh,plotflag)
%% Calculates waveform width and FWHM and automatically assigns
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% fwThresh = low cutoff for waveforms Full Width (arbitrary units)
% fwhmThresh = low cutoff for waveforms Full Width Half Maximum (arbitrary units)
% frThresh = low cutoff for waveforms Firing Rate (Hz)
% plotflag = binary of whether to plot the output
%
% Outputs:
% root = modified root object with root.info.uType putative unit type assignment
%       1 = Putative pyramidal
%       0 = Putative IN
%
% Created 3/19/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    fwThresh = 15
    fwhmThresh = 5
    frThresh = 100
    plotflag = 1
end

unitFRs  = [];
unitAmps = [];
unitFWHM = [];
unitFW   = [];
unitProm = [];
unitcls  = [];
unitprhp = []; % pre-hyperpolarization
unitpohp = []; % post-hyperpolarization

for i = 1:height(root.info)
    % unit = root.info.cluster_id == i;

    unitFRs  = [unitFRs; root.info.fr(i)];
    unitAmps = [unitAmps; root.info.Amplitude(i)];
    unitcls  = [unitcls; root.info.cluster_id(i)];

    % Extract waveform peak
    tmpwf = root.templateWF(i,:);
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

    unitprhp = [unitprhp; abs(mean(tmpwf) - firsthalf(lowind(1)))]; % pre-hyperpolarization
    unitpohp = [unitpohp; abs(mean(tmpwf) - lasthalf(lowind(2)))];  % post-hyperpolarization
    unitFWHM = [unitFWHM; tmpwidth(maxind)];
    unitFW   = [unitFW; (lowind(2) + length(firsthalf)) - lowind(1)];
    unitProm = [unitProm; tmpprom(maxind)];
end

%% WF Width cutoff

INWFs = unitFW <= fwThresh & unitFWHM <= fwhmThresh & unitFRs <= frThresh;
root.info.uType = ~INWFs;

% [goodcls(narrWFs), goodAmps(narrWFs), goodFW(narrWFs), goodFWHM(narrWFs)]

%% Plot for visual inspection

if plotflag
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

    figure;
    plot3(unitprhp(INWFs),unitFW(INWFs),unitFWHM(INWFs),'ro')
    hold on
    plot3(unitprhp(~INWFs),unitFW(~INWFs),unitFWHM(~INWFs),'bo')
    xlabel('pre-hyperpol')
    ylabel('FullWidth')
    zlabel('FWHM')

    goodINs = find(root.goodind & INWFs);
    goodPyrs = find(root.goodind & ~INWFs);
    nSubFigs = ceil(length(goodINs)/16);
    ct = 0;
    for i = 1:nSubFigs
        figure
        for j = 1:16
            try
                subplot(4,4,j)
                hold on
                plot(root.templateWF(goodINs(j+ct),:),'k')
                plot(root.templateWF(goodPyrs(j+ct),:),'b')
                plot([0 size(root.templateWF,2)],[0 0], 'k--')
                xlim([0 size(root.templateWF,2)])
                title(['Unit ' num2str(root.info.cluster_id(goodINs(j+ct)))])
                legend({['FR ' num2str(root.info.fr(goodINs(j+ct)))]},'Location','se','FontSize',8)
                set(gca,'FontSize',12,'FontName','Arial')
            catch
                continue
            end
        end
        ct = ct + 16;
    end
end

%% Optional plots for testing purposes
% figure;
% plot3(goodFRs,goodFW,goodFWHM,'ko')
% xlabel('mean FR')
% ylabel('goodFW')
% zlabel('FWHM')
% 
% figure;
% subplot(2,2,1)
% plot(goodProm./goodFWHM,goodFRs,'ko')
% ylabel('mean FR')
% xlabel('Prominence/FWHM ratio')
% xlim([0, max(goodProm./goodFWHM)])
% 
% subplot(2,2,2)
% plot(goodFWHM,goodFRs,'ko')
% ylabel('mean FR')
% xlabel('WF FWHM (AU)')
% xlim([0, max(goodFWHM)])
% 
% subplot(2,2,3)
% histogram(goodProm./goodFWHM,0:0.1:max(goodProm./goodFWHM))
% ylabel('Count')
% xlabel('Prominence/FWHM ratio')
% xlim([0, max(goodProm./goodFWHM)])
% 
% subplot(2,2,4)
% histogram(goodFWHM,0:0.5:max(goodFWHM))
% ylabel('Count')
% xlabel('WF FWHM (AU)')
% xlim([0, max(goodFWHM)])
% 
% figure;
% subplot(2,2,1)
% plot(goodAmps,goodFRs,'ko')
% ylabel('mean FR')
% xlabel('Amplitude(AU)')
% xlim([0, max(goodAmps)])
% 
% subplot(2,2,2)
% plot(goodFW,goodFRs,'ko')
% ylabel('mean FR')
% xlabel('WF Width (AU)')
% xlim([0, max(goodFW)])
% 
% subplot(2,2,3)
% histogram(goodAmps,0:0.5:max(goodAmps))
% ylabel('Count')
% xlabel('Amplitude (AU)')
% xlim([0, max(goodAmps)])
% 
% subplot(2,2,4)
% histogram(goodFW,0:1:max(goodFW))
% ylabel('Count')
% xlabel('WF Width (AU)')
% xlim([0, max(goodFW)])

%% PCA for waveform clustering -- doesn't seem to work consistently

% wfMetricsMat = [goodFRs,goodFW,goodFWHM,goodAmps,goodProm];
% 
% % [coefs,score,~,tsquared,explained] = pca(wfMetricsMat);
% 
% wfsMat = root.templateWF(root.goodind,:);
% [coefs,score,latent,tsquared,explained] = pca(wfsMat','centered',false);
% 
% figure; hold on;
% plot(explained,'k-o')
% plot(cumsum(explained),'k-o')
% 
% xlim([0 length(explained)])
% xlabel('PC#')
% ylabel('Explained variance')
% 
% figure; 
% plot3(score(:,1),score(:,2),score(:,3),'k.')
% xlabel('PC1')
% ylabel('PC2')
% zlabel('PC3')
% 
% hold on
% plot3(score(INWFs(root.goodind),1),score(INWFs(root.goodind),2),score(INWFs(root.goodind),3),'r.')
% title(replace(root.name,'_',' '))
% 
% figure; hold on
% plot(coefs(:,1),coefs(:,2),'.')
% plot(coefs(INWFs(root.goodind),1),coefs(INWFs(root.goodind),2),'r.')
% 
% % [goodcls, score(:,1), score(:,2), score(:,3)]
% 
% Z = linkage(coefs(:,1:3));
% figure
% dendrogram(Z)
