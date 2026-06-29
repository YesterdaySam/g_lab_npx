function [root,datStruc] = get_lyrBounds_ec(root,sess)
%% Assigns ec units based on peak on theta phase preference and/or opto tag
%
% Inputs:
% root = root object. Must have uPSDMax field
% minWidth = minimum layer width (100um) or use string to use peak-defined
% fields for all
%
% Outputs:
% root = updated root with root.lyrbounds
%
% Created 5/21/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------


tmpD = root.info.depth(root.goodind);

for i = 1:length(root.good)
    [datStruc.thetastats(i),datStruc.thFR(i,:)] = plot_thetaMod(root,root.good(i),1,2*pi/36,0);
    [datStruc.optobins,datStruc.optoMat(i,:,:)] = plot_frXopto(root,root.good(i),sess,0.002,0.5,0);
    datStruc.optoPkT(i) = get_firstPk(datStruc,i);
end
[datStruc.thAng, ~, datStruc.thSigs] = get_thAng(datStruc.thetastats);
datStruc.thSig = datStruc.thSigs < 0.05;
% shortLat = datStruc.optoPkT <= 0.01;

%% Assign putative layer - using theta angle

ec2inds = datStruc.thAng > -90 & datStruc.thAng < 90;
ec3inds = ~ec2inds; %-180:-90 and 90:180;

dBins = 0:50:max(root.info.depth);
lyrct = histcounts(tmpD,dBins); % All good by depth bin
lyrsg = histcounts(tmpD(datStruc.thSig),dBins); % Significantly theta by depth
lyrtr2 = histcounts(tmpD(ec2inds),dBins); % Theta tuning by depth
ec5stt = dBins(find(smooth(lyrsg./lyrct) <= 0.1,1)); % Estimate L5 start using cutoff of ratio of significant / all
ec3stt = dBins(find(smooth(lyrtr2./lyrct) <= 0.4,1)); % Estimate L3 start using cutoff of ratio of L2 theta / all

ec5inds = tmpD' > ec5stt & ~datStruc.thSig;
ec2inds = ec2inds & tmpD' < ec3stt + 300;
try
    ec3inds = (~ec2inds & ~ec5inds) | datStruc.optoPkT < 0.015;
catch
    ec3inds = ~ec2inds & ~ec5inds;
end

root.info.lyrID = zeros(height(root.info),1);
tmpInds = find(root.goodind);
notInds = find(~root.goodind);
ec2all  = root.info.depth(notInds) < ec3stt;
ec3all  = root.info.depth(notInds) > ec3stt;
ec5all  = root.info.depth(notInds) > ec5stt;
root.info.lyrID(tmpInds(ec2inds)) = 2;
root.info.lyrID(notInds(ec2all))  = 2;
root.info.lyrID(tmpInds(ec3inds)) = 3;
root.info.lyrID(notInds(ec3all))  = 3;
root.info.lyrID(tmpInds(ec5inds)) = 5;
root.info.lyrID(notInds(ec5all))  = 5;
root.lyrbounds = [ec3stt; ec5stt];

datStruc.lyrID = zeros(size(root.good));
datStruc.lyrID(ec5inds) = 5;
datStruc.lyrID(ec2inds) = 2;
datStruc.lyrID(ec3inds) = 3;

%% Assign putative layer - using opto
% 
% % moving average of opto pulse plus low latency units
% [sortDepth,sortInd] = sort(tmpD);
% sortOpto = datStruc.optoPkT(sortInd);
% meanOpPk = smoothdata(sortOpto,'movmean',10);
% ec3stt = sortDepth(find(meanOpPk < 0.025,1,'first'));
% ec3end = sortDepth(find(meanOpPk < 0.025,1,'last'));
% 
% % figure; hold on
% % plot(sortOpto)
% % plot(meanOpPk)
% % xlabel('unit (sorted by depth)')
% % ylabel('First opto peak post-pulse (s)')
% % legend('Raw','Moving Mean')
% % set(gca,'FontSize',12,'FontName','Arial')
% 
% ec5inds = tmpD > ec3end;
% ec2inds = tmpD < ec3stt;
% ec3inds = (tmpD >= ec3stt & tmpD <= ec3end) | datStruc.optoPkT' < 0.015;
% 
% datStruc.lyrID = zeros(size(root.good));
% datStruc.lyrID(ec5inds) = 5;
% datStruc.lyrID(ec2inds) = 2;
% datStruc.lyrID(ec3inds) = 3;

% opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time
% plot(datStruc.optoPkT(datStruc.lyrID == 5)*1000, tmpD(datStruc.lyrID == 5),'k*')
% plot(datStruc.optoPkT(datStruc.lyrID == 3)*1000, tmpD(datStruc.lyrID == 3),'c*')
% plot(datStruc.optoPkT(datStruc.lyrID == 2)*1000, tmpD(datStruc.lyrID == 2),'m*')
% % legend('Data','EC5','EC3','EC2')
% set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]); xlim([0 100]); xticks([0 100]); xticklabels([0 100])
% 
% thDepthFig = plot_datXdepth(root,datStruc,1,0,0,2); % theta
% plot(datStruc.thAng(datStruc.lyrID == 5)+180, tmpD(datStruc.lyrID == 5),'k*')
% plot(datStruc.thAng(datStruc.lyrID == 3)+180, tmpD(datStruc.lyrID == 3),'c*')
% plot(datStruc.thAng(datStruc.lyrID == 2)+180, tmpD(datStruc.lyrID == 2),'m*')
% % legend('Data','EC5','EC3','EC2')
% set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]);
% 
% saveas(opDepthFig,[sbase '_opXdepth.png'])
% saveas(thDepthFig,[sbase '_thXdepth.png'])

% thDepthFig = plot_datXdepth(root,datStruc,1,0,0,2); % theta
% plot(datStruc.thAng(datStruc.lyrID == 5)+180, tmpD(datStruc.lyrID == 5),'k*')
% plot(datStruc.thAng(datStruc.lyrID == 3)+180, tmpD(datStruc.lyrID == 3),'c*')
% plot(datStruc.thAng(datStruc.lyrID == 2)+180, tmpD(datStruc.lyrID == 2),'m*')
% xlim([0 360]); ylabel('Distance from tip (mm)'); xlabel('Theta angle')
% % legend('Data','EC5','EC3','EC2')
% % set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6]);

end