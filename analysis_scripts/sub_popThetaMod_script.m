%% Subiculum theta-modulation at population level

spath = 'D:\Data\Kelton\analyses\KW032\KW032_03122025_rec_D6_LMed1';
saveFlag = 1;

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)

nshanks = numel(unique(root.lfpinfo.lfpShank));

%% Break out LFP for old root files for easier subsequent loading

if saveFlag
    lfp.name     = root.name;
    lfp.lfp      = root.lfp;
    lfp.fs_lfp   = root.fs_lfp;
    lfp.lfpinfo  = root.lfpinfo;
    lfp.lfp_tsb  = root.lfp_tsb;
    lfp.bands    = root.bands;
    lfp.uPSDMax  = root.uPSDMax;
    lfp.uPSD     = root.uPSD;

    root.lfp     = root.lfp(root.uPSDMax(2,:),:);
    root.uPSD    = root.uPSD(:,root.uPSDMax(2,:));
    root.lfpinfo = root.lfpinfo(root.uPSDMax(2,:),:);
    root.uPSDMax = [repmat(1:4,size(root.bands,2),1)]; % Set uPSD to match updated LFP info

    for i = 1:nshanks
        thetaLFPSh = bandpass(root.lfp(root.uPSDMax(2,i),:),root.bands(1,:), root.fs_lfp);
        root.thEnv(i,:) = hilbert(thetaLFPSh);  % Operates along rows
    end

    saveRoot(root,spath)
    save([lfp.name '_lfp'],'lfp')
end

%% Collect all good, in-layer unit spikes relative to local or same-reference theta

thbnsz = 2*pi/36;
phase_edges = 0:thbnsz:2*pi;
lfp_use = [1 2 3 4];    % Set per shank or use same shank for global ref

thMap = [];    % [cluster_id, shankID, P(ripple particip.), Sig. Modulation]
thInfo = [];
ct = 1;

for j = 1:nshanks
    th_phase = angle(root.thEnv(lfp_use(j),:))+pi;
    % th_occ = histcounts(th_phase,phase_edges)/root.fs_lfp;

    lyrCCs = find(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1);
    for i = 1:length(lyrCCs)
        cc = root.info.cluster_id(lyrCCs(i));

        circ_stats = plot_thetaMod(root,cc,lfp_use(j),thbnsz,0);
        spkInds = root.tsb(root.cl == cc);

        lfpCounts = histcounts(spkInds,root.lfp_tsb);

        inds = [];

        % find index of all bins >i spikes and concatenate to existing indices
        for k = 0:max(lfpCounts)
            inds = cat(2,inds,find(lfpCounts>k));
        end
        % sort inds in case of repeats
        lfpInds = sort(inds);

        thMap(ct,:) = histcounts(th_phase(lfpInds),phase_edges);
        thInfo(ct,:) = [cc, j-1, root.info.uType(lyrCCs(i)), circ_stats.mrl, circ_stats.ang, circ_stats.p];
        ct = ct + 1;
    end
end

%% Observational statistics

disp(['Overall in-layer significant Theta-Modulated: ' num2str(sum(thInfo(:,6) < 0.05)) '/' num2str(size(thInfo,1))])

for i = 1:nshanks
    disp(['Shank ' num2str(i-1) ' significant ' num2str(sum(thInfo(thInfo(:,2) == i-1,6) < 0.05)) '/' num2str(numel(thInfo(thInfo(:,2) == i-1,6)))])
    muMRL(i) = mean(thInfo(thInfo(:,2) == i-1,4));
    disp(['      ' num2str(i-1) ' mean resultant length: ' num2str(muMRL(i))])
    muAng(i) = rad2deg(mean(thInfo(thInfo(:,2) == i-1,5)));
    disp(['      ' num2str(i-1) ' mean angle: ' num2str(muAng(i))])
end

%% Plot for each shank
cmapcool = cool(nshanks);
legCell = {};

thModFig = figure; hold on
xcoords = rad2deg(phase_edges);
xcoords = [xcoords, xcoords(2:end) + 360];

for j = 1:nshanks
    shInds = thInfo(:,2) == j-1;
    if sum(shInds) < 2
        continue
    end
    tmpMap = [thMap(shInds,:), thMap(shInds,:), thMap(shInds,1)];
    tmpMapNorm = mean(tmpMap)./sum(mean(tmpMap),'all');
    tmpSD = std(tmpMap)./sum(std(tmpMap),'all');
    upCI = tmpMapNorm + 1.96*tmpSD/sqrt(numel(shInds));
    dnCI = tmpMapNorm - 1.96*tmpSD/sqrt(numel(shInds));    
    % patch([xcoords fliplr(xcoords)],[dnCI fliplr(upCI)],cmapcool(j,:),'FaceAlpha',.1,'EdgeColor',cmapcool(j,:),'HandleVisibility','off')
    % plot(xcoords,tmpMapNorm,'Color',cmapcool(j,:))
    errorbar(xcoords,tmpMapNorm,[],tmpSD/sqrt(numel(shInds)),'Color',cmapcool(j,:))
    plot(muAng(j)+360,max(tmpMapNorm)+0.005,'*','Color',cmapcool(j,:),'HandleVisibility','off')
    legCell{j} = ['Shank ' num2str(j-1)];
end
plot(xcoords,-cos([phase_edges phase_edges(1:end-1)])*0.001+0.008,'k')
ylabel('Normalized Spike Probability'); ylims = ylim; ylim([ylims(1)-0.1*ylims(1), ylims(2)+0.1*ylims(2)]);
xlabel('Theta Phase'); xlim([0 720]); xticks(linspace(0, 720, 9));
legend(legCell,'Location','southeast')
set(gca,'FontSize',12,'FontName','Arial')

saveas(thModFig,[root.name '_ProbThSpike.png'])

%% Save for easier loading later

save([root.name '_ProbThSpike_data'], 'thMap', 'thInfo')








      