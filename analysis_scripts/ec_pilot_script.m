% EC
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05062025_rec_D1_RMed1';
spath = 'D:\Data\Kelton\analyses\KW043\KW043_05072025_rec_D2_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05082025_rec_D3_RLat2';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05092025_rec_D4_RMed2';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)
epochfile = dir("*_ec_pilot_dat.mat");
try load(epochfile.name); catch; disp('No existing epoched data file'); end

rwdShift = sess.valTrials(find(diff(sess.pos(sess.rwdind)) > 0.1,1)+1);   % Find lap of reward shift
nUnits = length(root.good);
saveFlag = 1;

dbnsz = 0.05;
histoBnsz = 5;
wlen = 150;
ripRef = root.ripRef;
r1pos = 0.9;    % 10 cm
binpos = 0.025:dbnsz:1.825;
tmpD = root.info.depth(root.goodind);

sbase = root.name; 

%% Plot example pre/post behavior

lickDIFig = plot_lickDiscrim(sess,[r1pos 1.8]*100,30,10);

if saveFlag
    saveas(lickDIFig,[sbase '_lickDI'], 'png')
end

%% Get real unit and epoch parameters

clear datStruc
datStruc.name   = root.name;
datStruc.trueSI = zeros(size(root.good));
datStruc.truePk = zeros(size(root.good));
datStruc.trueLc = zeros(size(root.good));

for i = 1:nUnits
    cc = root.good(i);
    lfpInd = root.info.shankID(root.info.cluster_id == cc)+1;   % Account for 0-indexing
    [datStruc.trueSI(i),~,datStruc.truePk(i),datStruc.trueLc(i),~,~,datStruc.posfr(i,:),datStruc.binedges] = get_SI(root,cc,sess,dbnsz);
    datStruc.frStandRun(i,:) = get_frStandVRun(root,cc,sess);
    [~,~,datStruc.trueVelMdl(i)] = plot_frXvel(root,cc,sess,2,0);
    [datStruc.thetastats(i),datStruc.thetafr(i,:)] = plot_thetaMod(root,cc,lfpInd,2*pi/36,0);
    % datStruc.swrfr(i,:) = plot_frXripple(root,cc,sess,ripRef,wlen,histoBnsz,0);
    [~,datStruc.rwdfr(i,:),datStruc.trueRI(i)] = plot_frXrwdtime(root,cc,sess,0.25,5,0);
    [datStruc.optobins,datStruc.optoMat(i,:,:)] = plot_frXopto(root,cc,sess,0.002,0.5,0);
end

datStruc = subEpochSI(sess,root,datStruc,10);

% datStruc.ripRate = size(root.ripStruc(ripRef).ripples,1) / (sum(not(sess.runInds)) / sess.samprate);    %Normalize based on standing periods
datStruc.binpos = datStruc.binedges(1:end-1)+0.5*dbnsz;

%% Find place cells and SPW-R modulation pre and post shift
% Methods: 
% Kitanishi et al. 2021: Exceed 99th percentile of SI shuffle
% Grienberger & Magee 2022:
%   Field: contiguous area with >= 20% peak FR
%   Induction: First lap with >3SD of in-field FR, above out-field noise
%              and activity in 2/5 subsequent laps also >3SD
%   Stability: Significant >3SD activity in 30% of post-induction laps

% Add method for jittering opto tagging

tic
nShufs = 250;

datStruc.shufSI = zeros(nUnits,nShufs);
datStruc.shufSPWR = zeros(nShufs,nUnits,length(-wlen:histoBnsz:wlen)-1);

for j = 1:nShufs
    [shiftroot,shiftsess,shiftInd] = shiftTrain(root,sess);

    if mod(j,50) == 0
        disp(['Shuffle # ' num2str(j)])
        toc 
    end

    for i = 1:nUnits
        cc = root.good(i);

        [datStruc.shufSI(i,j)] = get_SI(shiftroot,cc,shiftsess);

        % frstHalf.shufSPWR(j,i,:) = plot_frXripple(shiftFrst,cc,shiftSessFrst,ripRef,wlen,histoBnsz,0);

        % [~,~,frstHalf.shufRI(i,j)] = plot_frXrwdtime(shiftFrst,cc,shiftSessFrst,0.25,5,0);
        % [~,~,lastHalf.shufRI(i,j)] = plot_frXrwdtime(shiftLast,cc,shiftSessLast,0.25,5,0);
    end
end

toc

datStruc.sig = sum(datStruc.shufSI > datStruc.trueSI,2) / nShufs;

disp(['Sig SI p <= 0.05 ' num2str(sum(datStruc.sig <= 0.05)) ' of ' num2str(nUnits) ' units'])

%% Save for later

if saveFlag
    save([root.name '_ec_pilot_dat'],'datStruc')
end

%% Future analyses
% Plot spatial frequency of each cell for 1D grid cells

for i = 1:length(root.good)
    tmp = get_firstPk(datStruc,i);
    datStruc.optoPkT(i) = tmp;
end

thDepthFig = plot_datXdepth(root,datStruc,1,0,0,2); % theta
siDepthFig = plot_datXdepth(root,datStruc,1,0,0,3); % spatial info
opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time

if saveFlag
    saveas(thDepthFig,[sbase '_ec_pilot_thXdepth.png'])
    saveas(siDepthFig,[sbase '_ec_pilot_siXdepth.png'])
    saveas(opDepthFig,[sbase '_ec_pilot_opXdepth.png'])
end

%% Assign putative layer

for i = 1:length(root.good)
    datStruc.thAng(i) = rad2deg(datStruc.thetastats(i).ang);
    datStruc.thSig(i) = datStruc.thetastats(i).p < 0.05;
end

shortLat = datStruc.optoPkT <= 0.01;

figure; hold on
plot(datStruc.thAng+180,datStruc.optoPkT,'ko')
plot(datStruc.thAng(shortLat)+180,datStruc.optoPkT(shortLat),'r*')

% PCA attempt
% clusterDat = [normalize(datStruc.thAng'+180,'range'),normalize(datStruc.optoPkT','range'),normalize(tmpD,'range')];
% 
% [coeff,score,latent,tsquare,explained] = pca(clusterDat);
% 
% figure; 
% plot3(score(:,1),score(:,2),score(:,3),'k.')

% moving average of opto pulse plus low latency units
[sortDepth,sortInd] = sort(tmpD);
sortOpto = datStruc.optoPkT(sortInd);
meanOpPk = smoothdata(sortOpto,'movmean',10);
ec3stt = sortDepth(find(meanOpPk < 0.025,1,'first'));
ec3end = sortDepth(find(meanOpPk < 0.025,1,'last'));

figure; hold on
plot(sortOpto)
plot(meanOpPk)
xlabel('unit (sorted by depth)')
ylabel('First opto peak post-pulse (s)')
legend('Moving Mean','Raw')
set(gca,'FontSize',12,'FontName','Arial')

ec5inds = tmpD > ec3end;
ec2inds = tmpD < ec3stt;
ec3inds = (tmpD >= ec3stt & tmpD <= ec3end) | datStruc.optoPkT' < 0.015;

datStruc.lyrID = zeros(size(datStruc.trueSI));
datStruc.lyrID(ec5inds) = 5;
datStruc.lyrID(ec2inds) = 2;
datStruc.lyrID(ec3inds) = 3;

opDepthFig = plot_datXdepth(root,datStruc,1,0,0,4); % peri-opto peak time

plot(datStruc.optoPkT(datStruc.lyrID == 5)*1000, tmpD(datStruc.lyrID == 5),'k*')
plot(datStruc.optoPkT(datStruc.lyrID == 3)*1000, tmpD(datStruc.lyrID == 3),'c*')
plot(datStruc.optoPkT(datStruc.lyrID == 2)*1000, tmpD(datStruc.lyrID == 2),'m*')

%% Descriptive statistics and graphs
% root.good(sigFrst <= 0.05 & sigLast <= 0.05 & root.info.fr(root.goodind)' > 0.1 & root.info.lyrID(root.goodind)' == 1)

% lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = hiFRUnits & root.info.uType(root.goodind);   % hi FR and not IN
siUnits = useUnits & datStruc.sig < 0.05;

nUseUnits = sum(useUnits);

xcoords = ones(nUseUnits,1);
shiftbins = [-fliplr(datStruc.binedges) datStruc.binedges(2:end)];

% Firing rate stand vs run
[~,ps.standFR_firstlast] = ttest(datStruc.frStandRun(useUnits,1),lastHalf.frStandRun(useUnits,1));
[~,ps.runFR_firstlast] = ttest(datStruc.frStandRun(useUnits,2),lastHalf.frStandRun(useUnits,2));

% Spatial Information
[~,ps.SI_firstlast] = ttest(datStruc.trueSI(siUnits),lastHalf.trueSI(siUnits));

lcBothSIFig = plotBar2(datStruc.trueSI(siUnits), lastHalf.trueSI(siUnits));
ylabel('Spatial Information (Bits/spike)')

if saveFlag
    sbase = root.name;
    saveas(frStndRunFig,[sbase '_RwdShift_FRStandRun.png'])
    saveas(lcBothSIFig,[sbase '_RwdShift_SI.png'])
end

%% Peak distribution
datStruc.rwdBin = find(datStruc.binpos > r1pos,1);

ecFieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(siUnits)),datStruc.binedges);

[spWflFig] = plot_unitsXpos(root,sess,root.good(siUnits));
plot([datStruc.rwdBin datStruc.rwdBin],[0 length(siUnits)+1],'r--','LineWidth',2)
% plot(normalize(ecFieldDistro,'range').*sum(siUnits)./5,'k','LineWidth',1)
title(replace(root.name,"_"," "))

peakDistroFig = plotDistroHisto(ecFieldDistro,binpos,r1pos);
set(gcf,'units','normalized','position',[0.4 0.35 0.3 0.14])
xlabel('Track Position (cm)')

if saveFlag
    saveas(spWflFig,[sbase '_ec_pilot_waterfall.png'])
    saveas(peakDistroFig,[sbase '_ec_pilot_PeakDistro.png'])
end

%% Functions

function [firstPkT] = get_firstPk(datStruc,unit)

tmpOpto = squeeze(datStruc.optoMat(unit,:,:));

baseInds = datStruc.optobins < 0;
counts = sum(tmpOpto);
baseMu = mean(counts(baseInds));
baseSD = std(counts(baseInds));

countZ = (counts - baseMu) ./ baseSD;
smoothZ = smoothdata(countZ,'gaussian',10);

[tmppks, tmplocs] = findpeaks(smoothZ);

% figure; hold on
% plot(datStruc.optobins,countZ)
% plot(datStruc.optobins,smoothZ)
% plot(datStruc.optobins(tmplocs),tmppks,'v')

% first peak >2 Z score
postPulsePks = datStruc.optobins(tmplocs) > 0;
hiPks = tmppks > 2;

firstPkT = datStruc.optobins(tmplocs(find(postPulsePks & hiPks, 1)));

if isempty(firstPkT)
    firstPkT = NaN;
end

end

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
