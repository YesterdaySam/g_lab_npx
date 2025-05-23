%% Script to generate a crude peri-opto plot by units
optoTS = sess.ts(sess.optoind);
ctlTS = sess.ts(root.tsb(root.cl == 57)); %EC4= 57; EC3 = 125
% dShift = ctlTS - optoTS(361:540);     %For EC D2 bank 3 / EC3
dShift = ctlTS - optoTS(548:720);     %For EC D2 bank 4 / EC4
disp(num2str(mode(dShift)))

optoShiftTs = optoTS + mode(dShift);

nUnits = length(root.good);

bnsz = 0.002;
bnrg = 0.1;

%%
for i = 1:nUnits
    % [bnctrs,tmpbnfr] = plot_frXopto_shift(root,optoShiftTs,root.good(i),sess,0.002,0.1,0);
    [bnctrs,tmpbnfr] = plot_frXopto(root,root.good(i),sess,0.002,0.1,0);
    frMapRaw(i,:) = mean(tmpbnfr,1,'omitnan');
end
nBins = size(frMapRaw,2);


unitMax = max(frMapRaw,[],2);
frMapNorm = frMapRaw ./ repmat(unitMax,[1, nBins]);
% lininds = find(frMapNorm == 1); %Returns linear indexing of MxN matrix
% [is,js] = ind2sub(size(frMapNorm), lininds);
for i = 1:nUnits
    maxBin(i) = find(frMapNorm(i,:) == 1,1);
end
[~,sortInd] = sort(maxBin);

frMapSort = frMapNorm(sortInd,:);

fhandle = figure; hold on;
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.5])
imagesc(frMapSort,[prctile(frMapSort,1,'all'), prctile(frMapSort,98,'all')]);
plot([50 50],[0 nUnits+1],'w--','LineWidth',2)
colormap("copper")
cbar = colorbar; clim([0 0.98]);
xticks(0:10:nBins+1)
xticklabels((-0.1:0.02:0.1)*1000)
xlim([0.15 nBins+0.85]); ylim([0.15 nUnits+0.85])
xlabel('Time relative to opto pulse (ms)');
ylabel('Unit #'); ylabel(cbar,'FR (Normalized)','FontSize',12,'Rotation',90)
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

%% Store vars by session
frMapSort = [frMapSort_EC3; frMapSort_EC4];

nUnits = size(frMapSort,1);

%% Split units by opto response

respY = [];
respN = [];

for i = 1:nUnits
    meanFR = mean(frMapSort(i,:));
    stdevFR = std(frMapSort(i,:));
    maxBin = find(frMapSort(i,:) == 1,1);
    if 1 > meanFR + 5*stdevFR & maxBin > nBins / 2 + 1
        respY = [respY; frMapSort(i,:)];
    else
        respN = [respN; frMapSort(i,:)];
    end
end

%% Plot 

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.2])
% for i = 1:size(respY,1)
%     plot(bnctrs,smooth(respY(i,:),5),'Color',[0.5,0.5,0.5]);
% end
plot(bnctrs*1000, respY, 'Color', [0.5,0.5,0.5])
plot(bnctrs*1000, smooth(mean(respY,1),3), 'k', 'LineWidth',2)
xlabel('Time relative to opto pulse (ms)'); ylabel('FR (Normalized)')
set(gca,'FontSize',12,'FontName','Arial')

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.2])
plot(bnctrs*1000, respN, 'Color', [0.5,0.5,0.5])
plot(bnctrs*1000, smooth(mean(respN,1),3), 'k', 'LineWidth',2)
xlabel('Time relative to opto pulse (ms)'); ylabel('FR (Normalized)')
set(gca,'FontSize',12,'FontName','Arial')
