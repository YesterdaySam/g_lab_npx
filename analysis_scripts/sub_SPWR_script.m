%% 
% Key Methods
% Pointwise/Global confidence bands Fujisawa et al., 2008; 10.1038/nn.2134
% Ripple-Spike CCGs Oliva et al., 2016; 10.1016/j.neuron.2016.08.008
% Spike jitter method Amarasingham et al., 2011; 10.1152/jn.00633.2011

spath = 'D:\Data\Kelton\analyses\KW040\KW040_04292025_rec_D4_LLat1';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)

rwdShift = sess.valTrials(find(diff(sess.pos(sess.rwdind)) > 0.1)+1);   % Find lap of reward shift

%% Make summary plots

tmpRipDur = plot_ripDurHisto(root,sess,5,150);
tmpRipPwr = plot_ripPwrHisto(root);
[tmpRipIRI, tmpRipIRIZoom] = plot_ripIRIHisto(root,sess);
tmpRipTime = plot_ripTimeHisto(root,sess);
if ~isempty(rwdShift)
    plot([sess.ts(sess.lapstt(rwdShift)),sess.ts(sess.lapstt(rwdShift))],[0 25],'k--','LineWidth',2,'HandleVisibility','off');
end
    
if saveFlag
    saveas(tmpRipDur,[root.name '_rippleDur.png'])
    saveas(tmpRipPwr,[root.name '_ripplePwr.png'])
    saveas(tmpRipIRI,[root.name '_rippleIRI.png'])
    saveas(tmpRipIRIZoom,[root.name '_rippleIRI_zoom.png'])
    saveas(tmpRipTime,[root.name '_rippleTiming.png'])
end

%% Batch plot units relative to ripples
wlen = 150;
histoBnsz = 5;
ripref = 2; % Index of ripStruc to use as reference ripples
if ~isempty(rwdShift)
    postShiftInd = find(sess.ts(root.ripStruc(ripref).ripples(:,3)) > sess.ts(sess.lapstt(rwdShift)),1);
end

mkdir("ripplePlots_good")
cd('ripplePlots_good')

for i = 1:length(root.good)
    cc = root.good(i);
    [~,~,tmpfig] = plot_frXripple(root,cc,sess,ripref,wlen,histoBnsz);
    tmpRipPart = get_RipParticipation(root,cc,sess,ripref,wlen);

    if ~isempty(rwdShift)
        yyaxis left
        plot([-wlen wlen],[postShiftInd postShiftInd],'k--')    
    end

    legend(['P(Participation) ' num2str(tmpRipPart)])
    tmpind = find(root.info.cluster_id == cc);
    saveas(tmpfig, ['Unit' num2str(cc) '_Shank' num2str(root.info.shankID(tmpind)), '_Depth', num2str(root.info.depth(tmpind)) '_ripRaster'], 'png')
    close(tmpfig)
end
cd('..')

%% Test 1 unit against 95% Global Error band shuffle (circular or jitter)
cc = 171;
nShufs = 500;
jitlen = round(150/1000/(1/sess.samprate)); % -100:100 msec jitter interval
histoBnsz = 5;
wlen = 150;

[orig_count,binedges] = plot_frXripple(root,cc,sess,ripref,wlen,histoBnsz,1);
% tmpRipMod = get_RipMod(root,cc,sess,1000);
tmpRipPart = get_RipParticipation(root,cc,sess,ripref,wlen);
legend(['P(Participation) ' num2str(tmpRipPart)])
jit_count = zeros(nShufs,length(orig_count));

for i = 1:nShufs
    % Circular shuffle method
    shiftRoot = shiftTrain(root,sess);
    jit_count(i,:) = plot_frXripple(shiftRoot,cc,sess,ripref,wlen,histoBnsz,0);

    % % Jitter shuffle method
    % jits = jitter_spks(root.tsb(root.cl == cc), jitlen);   
    % % jits = jitter_spks(cc, jitlen);   
    % jit_count(i,:) = plot_frXripple(root,jits,sess,wlen,histoBnsz,0);
end

[~,tmpfig] = get_confband(jit_count,orig_count,1,binedges(1:end-1)*1000,histoBnsz/2);
xlabel('Time to ripple peak (sec)')

if saveFlag
    tblind = find(root.info.cluster_id == cc);
    sbase = ['unit_' num2str(cc) '_shank_' num2str(root.info.shankID(tblind)) '_rwdshift_'];
    saveas(tmpfig,[sbase 'spwrMod'], 'png')
end

%% Test all good units in layer for modulation
nShufs = 500;
jitlen = round(150/1000/(1/sess.samprate)); % -100:100 msec jitter interval
histoBnsz = 5;
wlen = 150;

ripMod = [];    % [cluster_id, shankID, P(ripple particip.), Sig. Modulation]
ct = 1;

for j = 1:nshanks
    root.ripples = ripStruc(j).ripples;
    lyrCCs = find(root.goodind & root.info.shankID == j-1 & root.info.lyrID == 1);
    for i = 1:length(lyrCCs)
        cc = root.info.cluster_id(lyrCCs(i));
        [orig_count,binedges] = plot_frXripple(root,cc,sess,ripref,wlen,histoBnsz,0);
        ripMod(ct,1:3) = [cc, j-1 root.info.uType(lyrCCs(i))];
        ripMod(ct,4) = get_RipParticipation(root,cc,sess,ripref,wlen);

        jit_count = zeros(nShufs,length(orig_count));
        for k = 1:nShufs
            % Circular shuffle method
            shiftRoot = shiftTrain(root,sess);
            jit_count(k,:) = plot_frXripple(shiftRoot,cc,sess,ripref,wlen,histoBnsz,0);

            % % Jitter shuffle method
            % jits = jitter_spks(root.tsb(root.cl == cc), jitlen);
            % % jits = jitter_spks(cc, jitlen);
            % jit_count(i,:) = plot_frXripple(root,jits,sess,wlen,histoBnsz,0);
        end

        p_final = get_confband(jit_count,orig_count);

        ripMod(ct,5) = sum(p_final);
        ripMod(ct,6) = sum(p_final) > 1;    % If more than 1 point falls outside global significance
        ct = ct + 1;
    end
    disp(['Finished shuffling shank ' num2str(j)])
end

root.ripMod = ripMod;

%% Observational statistics

disp(['Overall in-layer significant ' num2str(sum(ripMod(:,6))) '/' num2str(size(ripMod,1))])

for i = 1:nshanks
    disp(['Shank ' num2str(i-1) ' significant ' num2str(sum(ripMod(ripMod(:,2) == i-1,6))) '/' num2str(numel(ripMod(ripMod(:,2) == i-1,6)))])
    disp(['      ' num2str(i-1) ' mean participation: ' num2str(mean(ripMod(ripMod(:,2) == i-1,4)))])
end

%% Waterfall
histoBnsz = 1;
dummy_counts = plot_frXripple(root,1,sess,ripref,wlen,histoBnsz,0);
gaussSigma = 5;
sz = 32;
gaussX = linspace(-sz/2, sz/2, sz);
gaussFilt = exp(-gaussX .^2 / (2*gaussSigma ^2));
gaussFilt = gaussFilt / sum(gaussFilt); % Normalize

ct = 1;
ripZMap = [];

for i = 1:size(ripMod,1)
    if ripMod(i,6) == 1
        cc = ripMod(i,1);
        [tmpPeriSpikes,bins] = plot_frXripple(root,cc,sess,ripref,wlen,histoBnsz,0);
        muRipSpks = mean(tmpPeriSpikes);
        sdRipSpks = std(tmpPeriSpikes);
        zRipSpks = (tmpPeriSpikes - muRipSpks) ./ sdRipSpks;
        ripZMap(ct,:) = conv(zRipSpks,gaussFilt,'same');
    else
        ripZMap(ct,:) = zeros(1,length(dummy_counts));
    end
    ct = ct + 1;
end

sigMod = ripMod(:,6) == 1;
sigInfo = ripMod(sigMod,:);
sigMap = ripZMap(sigMod,:);

% Sort map by shank
for i = 1:nshanks
    inds = sigInfo(:,2) == i-1;
    tmpHistos = sigMap(inds,:);
    tmpInfo = sigInfo(inds,:);

    [maxZ, maxInds] = max(tmpHistos,[],2);    % Account for +/- modulated units
    [minZ, minInds] = min(tmpHistos,[],2);
    negInds = abs(minZ) > maxZ;
    maxInds(negInds) = minInds(negInds);

    [~,sortInds] = sort(maxInds);   % Sort by timing of peak, not peak height
    
    sigMap(inds,:) = tmpHistos(sortInds,:);
    sigInfo(inds,:) = tmpInfo(sortInds,:);
end

shankChgInds = find(diff(ripMod(sigMod,2)));

sigRipFig = figure;
imagesc(sigMap)
hold on
for i = 1:nshanks - 1
    plot([1 size(sigMap,2)],[shankChgInds(i)+.5 shankChgInds(i)+.5],'r--','LineWidth',2)
end
cbar = colorbar;
set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
ylabel('Unit #'); ylabel(cbar,'FR (Z-score)','FontSize',12,'Rotation',90)
xticks(linspace(0,size(sigMap,2),7))
xticklabels(linspace(-wlen,wlen,7))
xlabel('Time to Ripple Peak (msec)');

saveas(sigRipFig,[root.name '_sigRipModUnits.png'])

%% Save for easier loading later

save([root.name '_sigRipModUnit_data'], 'sigMap', 'sigInfo')


%% testing new shiftTrain.m

unit = 329;
origctall = numel(root.tsb(root.cl == unit));
origctrun = sum(sess.runInds(root.tsb(root.cl == unit)));

shufall = [];
shufrun = [];
for i = 1:100
    [tmproot,tmpsess] = shiftTrain(root,sess);
    shufall(i) = numel(root.tsb(root.cl == unit));
    shufrun(i) = sum(sess.runInds(root.tsb(root.cl == unit)));
end










