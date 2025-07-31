% EC
spath = 'D:\Data\Kelton\analyses\KW043\KW043_05062025_rec_D1_RMed1';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05072025_rec_D2_RLat1';
% spath = 'D:\Data\Kelton\analyses\KW043\KW043_05082025_rec_D3_RLat2';

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

sbase = root.name; 

%% Plot example pre/post behavior

lickDIFig = plot_lickDiscrim(sess,[r1pos 1.8]*100,30,10);

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

%% Descriptive statistics and graphs
% root.good(sigFrst <= 0.05 & sigLast <= 0.05 & root.info.fr(root.goodind)' > 0.1 & root.info.lyrID(root.goodind)' == 1)

% lyrUnits = root.info.lyrID(root.goodind) == 1;
hiFRUnits = root.info.fr(root.goodind) > 0.1;
useUnits = hiFRUnits & root.info.uType(root.goodind);
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
ecFieldDistro = histcounts(datStruc.binpos(datStruc.trueLc(siUnits)),datStruc.binedges);

peakDistroFig = plotDistroHisto(ecFieldDistro,binpos,r1pos);
xlabel('Track Position (cm)')

% % Align field distribution around reward
% trackLen = max(datStruc.binedges);
% shiftR1 = trackLen/2 - r1pos;
% fieldAlignR1 = mod(datStruc.binpos(datStruc.trueLc(siUnits))+shiftR1,trackLen)+0.5*dbnsz;    % Account for mod() creating 0's
% shiftR2 = trackLen/2 - r2pos;
% fieldAlignR2 = mod(lastHalf.binpos(lastHalf.trueLc(siUnits))+shiftR2,trackLen)+0.5*dbnsz;
% 
% fieldDistroFrst_Align = histcounts(fieldAlignR1,datStruc.binedges);
% fieldDistroLast_Align = histcounts(fieldAlignR2,lastHalf.binedges);
% 
% rzDistroFig = plotDistroHisto(fieldDistroFrst_Align,fieldDistroLast_Align,binpos-trackLen/2,0);
% xlabel('Distance to RZ (cm)')
% 
% Compare shift in PF peak relative to reward per unit
% deltaField_RZAlign = fieldAlignR2 - fieldAlignR1;
% circAlignNeg = deltaField_RZAlign < -trackLen/2;
% circAlignPos = deltaField_RZAlign > trackLen/2;
% deltaField_RZAlign(circAlignNeg) = deltaField_RZAlign(circAlignNeg) + trackLen;     % When new field forward-shifts
% deltaField_RZAlign(circAlignPos) = deltaField_RZAlign(circAlignPos) - trackLen;     % When new field back-shifts
% 
% deltaField_Distro = histcounts(deltaField_RZAlign,shiftbins);
% 
% fieldShiftRZFig = figure; hold on
% bar((shiftbins(1:end-1)+0.5*dbnsz)*100,deltaField_Distro/sum(deltaField_Distro),'FaceColor',[0.25 0.15 1])
% plot([0 0]*100,[0 0.15],'k--')
% xlabel('\Delta RZ-aligned Novel - Familiar (cm)')
% ylabel('Probability')
% xlim([-trackLen/2-dbnsz trackLen/2+dbnsz]*100)
% set(gca,'FontSize',12,'FontName','Arial')

if saveFlag
    saveas(peakDistroFig,[sbase '_ec_pilot_PeakDistro.png'])
    % saveas(rzDistroFig,[sbase '_RwdShift_PeakDistroRZ.png'])
    % saveas(fieldShiftRZFig,[sbase '_RwdShift_PeakShiftRZ.png'])
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
ylabel('Probability of field peak')
set(gca,'FontSize',12,'FontName','Arial')
end
