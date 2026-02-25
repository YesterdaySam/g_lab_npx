% DP2 progress report script

% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05062025_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05072025_rec_D2_RLat1';
spath = 'D:\Data\Kelton\analyses\KW043\KW043_05082025_rec_D3_RLat2';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05092025_rec_D4_RMed2';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_ec_pilot_dat.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

cd('D:\Data\Kelton\analyses\group_analyses\dp2_progreport')

%% Example opto tagging

ec2grp = [28, 85, 172, 131, 124];
ec3grp = [300, 384, 446, 451, 459];
ec5grp = [520, 546, 612, 665, 675];

ecgrp = ec3grp;

clear("normOptoFR")

for i = 1:length(ecgrp)
    tmpind = find(root.good == ecgrp(i));
    normOptoFR(i,:) = normalize(smoothdata(sum(squeeze(datStruc.optoMat(tmpind,:,:)),1),'gaussian',10),'norm');
end

ecOptoExFig = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.2])

for i = 1:length(ecgrp)
    plot(datStruc.optobins*1000, normOptoFR(i,:), 'Color', [.5 .5 .5])
end
ylabel('P(Spike)'); xlim([-250 250])
xlabel('Time to opto pulse (ms)')
plot([0 0], [0 0.3], 'k--')
set(gca,'FontSize',12,'FontName','Arial')

saveas(ecOptoExFig,'dp2_exOpto_ec3units','svg');

%% Modified opto x depth fig
% comment out extra overlay axis for shank in plot_datXdepth.m

for i = 1:length(ecgrp)
    ec2inds(i) = find(root.good == ec2grp(i));
    ec3inds(i) = find(root.good == ec3grp(i));
    ec5inds(i) = find(root.good == ec5grp(i));
end
tmpD = root.info.depth(root.goodind);

opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time
set(gcf,'units','normalized','position',[0.4 0.2 0.15 0.6])
xlim([0 250]); xticks([0 100 200]); xticklabels([0 100 200])
xlabel('Time from opto pulse (ms)'); hleg = legend(''); set(hleg,'visible','off')
plot(datStruc.optoPkT(ec2inds)*1000, tmpD(ec2inds),'m*')
plot(datStruc.optoPkT(ec3inds)*1000, tmpD(ec3inds),'c*')
plot(datStruc.optoPkT(ec5inds)*1000, tmpD(ec5inds),'k*')

saveas(opDepthFig,'dp2_depthXopto_kw043_D3','svg');

%% Waterfalls

% [fhandle,tmpMap,sortPre] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 3 & siUnits,:),datStruc.binedges);
ec3FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 3 & siUnits)),datStruc.binedges);
peak3DistroFig = plotDistroHisto(ec3FieldDistro,binpos,r1pos);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')

% [fhandle,tmpMap,sortPre] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 2 & siUnits,:),datStruc.binedges);
ec2FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 2 & siUnits)),datStruc.binedges);
peak2DistroFig = plotDistroHisto(ec2FieldDistro,binpos,r1pos);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')

% [fhandle,tmpMap,sortPre] = plot_unitWaterfall(datStruc.posfr(datStruc.lyrID == 5 & siUnits,:),datStruc.binedges);
ec5FieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(datStruc.lyrID == 5 & siUnits)),datStruc.binedges);
peak5DistroFig = plotDistroHisto(ec5FieldDistro,binpos,r1pos);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')




%%
function [fhandle] = plotDistroHisto(distro1,binpos,rzPos)
vColors2 = [0.5 0.5 1; 0.75 0.75 1];

fhandle = figure; hold on
plot(binpos*100,distro1./sum(distro1),'Color',vColors2(1,:),'LineWidth',2);
plot([rzPos rzPos]*100,[0 0.15],'k--')
ylim([0 max(distro1./sum(distro1),[],'all')])
xlim([binpos(1)*100-1 binpos(end)*100+1])
ylabel('P(field peak)')
set(gca,'FontSize',12,'FontName','Arial')
end


