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

%% compare shifted spike trains for pfs
cc = 28;

[tmpedges, tmpfr] = plot_frXpos(root,cc,sess,0.05,0.04,0);
[tmpSI] = get_SI(root,cc,sess);

for i = 1:1000
    shiftroot = shiftTrain(root,sess);
    [~,shufFR(i,:)] = plot_frXpos(shiftroot,cc,sess,0.05,0.04,0);
    [shufSI(i)] = get_SI(shiftroot,cc,sess);
    
end

figure; hold on
plot(tmpedges(1:end-1)*100,tmpfr,'b')
plot(tmpedges(1:end-1)*100,mean(shufFR),'k')
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
