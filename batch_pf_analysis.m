% root = root4;
% sess = sess4;

SI = zeros(1,length(root.good));
uFR = zeros(1,length(root.good));
pkFR = zeros(1,length(root.good));

nUnits = length(root.good);

for i = 1:nUnits
    cc = root.good(i);

    % [SI(i),uFR(i),pkFR(i)] = get_SI(root,cc,sess);
    [SI(i),pkFR(i),uFR(i)] = get_PF(root,cc,sess);
end

%% Find good PFs and Plot waterfalls of place fields
sithresh = 0.2;
frthresh = 2;
uThresh  = 10;

hiSI = SI > sithresh;
hiFR = pkFR > frthresh;
lowuFR = uFR < uThresh;

goodPFs = root.good(hiSI & hiFR & lowuFR);

tmpfig = plot_unitsXpos(root,sess,goodPFs);
% saveas(tmpfig,'goodSIFR_waterfall','png');

%% Find depth and shankID of good fields

for i = 1:numel(goodPFs)
    [goodDepth(i,1), goodDepth(i,2)] = getDepthByID(root,goodPFs(i));
end

%% compare shifted spike trains for pfs
cc = 207;
nShufs = 1000;

[tmpSI,tmppkFR,tmpuFR,~,tmpbinFR,tmpedges] = get_PF(root,cc,sess);
shufSI = zeros(nShufs,1);
shufuFR = zeros(nShufs,1);
shufpkFR = zeros(nShufs,1);
shufbinFR = zeros(nShufs,length(tmpedges)-1);

for i = 1:nShufs
    shiftroot = shiftTrain(root,sess);
    % [shufSI(i),shufpkFR(i),shufuFR(i),~,shufbinFR(i,:)] = get_PF(shiftroot,cc,sess); 
    [shufSI(i),shufuFR(i),shufpkFR(i),~,~,shufbinFR(i,:)] = get_SI(shiftroot,cc,sess);
end

figure; hold on
plot(tmpedges(1:end-1)*100,tmpbinFR,'b')
plot(tmpedges(1:end-1)*100,mean(shufbinFR),'k')
xlabel('Position (cm)')
ylabel('Firing Rate')
title(['Cell ' num2str(cc)])
legend({"True","Avg Shuffle"})
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

figure; hold on
histogram(shufSI,30)
if tmpSI > prctile(shufSI,99)
    plot([tmpSI, tmpSI], [0 100], 'r--')
else
    plot([tmpSI, tmpSI], [0 100], 'k--')
end
xlabel('Spatial Info.')
ylabel('Count')
title(['Cell ' num2str(cc)])
legend({"Shuffle","True"})
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')

shufTest = tmpSI > prctile(shufSI,99)

%% Batch shuffle all good units
nShufs = 1000;

shufTest = zeros(numel(root.good),1);

for i = 1:numel(root.good)
    cc = root.good(i);

    tmpSI = get_SI(root,cc,sess);
    shufSI = zeros(nShufs,1);

    for j = 1:nShufs
        shiftroot = shiftTrain(root,sess);
        [shufSI(j)] = get_SI(shiftroot,cc,sess);
    end

    shufTest(i) = tmpSI > prctile(shufSI,99);
end

figure; imagesc(shufTest)