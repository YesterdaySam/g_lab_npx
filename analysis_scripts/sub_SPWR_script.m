%% 
% Key Methods
% Pointwise/Global confidence bands Fujisawa et al., 2008; 10.1038/nn.2134
% Ripple-Spike CCGs Oliva et al., 2016; 10.1016/j.neuron.2016.08.008
% Spike jitter method Amarasingham et al., 2011; 10.1152/jn.00633.2011

spath = 'D:\Data\Kelton\analyses\KW022\KW022_12162024_rec_D1_RLat1';

cd(spath)
rootfile = dir("*_root.mat");
load(rootfile.name)
sessfile = dir("*_session.mat");
load(sessfile.name)

%% Get ripples at s. pyr. layer as determined by ripple band peak

nshanks = numel(unique(root.lfpinfo.lfpShank));

chans = root.uPSDMax(2,:);
nChans = length(chans);
cmapcool = cool(nChans);

ct = 1;
catRips = [];
for chan = chans
    ripStruc(ct).ripples = get_ripples(root,chan,sess,4,5,[15 250]);
    % catRips = [catRips; ripStruc(ct).ripples(:,2), -1+ct+zeros(size(ripStruc(ct).ripples,1),1)];
    ct = ct+1;
end

%% Batch plot units relative to ripples
wlen = 150;
histoBnsz = 5;

mkdir("ripplePlots_good")
cd('ripplePlots_good')

for i = 1:length(root.good)
    cc = root.good(i);
    [~,~,tmpfig] = plot_frXripple(root,cc,sess,wlen,histoBnsz);
    tmpRipPart = get_RipParticipation(root,cc,sess,wlen);
    legend(['P(Participation) ' num2str(tmpRipPart)])
    tmpind = find(root.info.cluster_id == cc);
    saveas(tmpfig, ['Unit' num2str(cc) '_Shank' num2str(root.info.shankID(tmpind)), '_Depth', num2str(root.info.depth(tmpind)) '_ripRaster'], 'png')
    close(tmpfig)
end
cd('..')

%% Test 1 unit against 95% Global Error band shuffle (circular or jitter)
cc = 245;
nShufs = 500;
jitlen = round(150/1000/(1/sess.samprate)); % -100:100 msec jitter interval
histoBnsz = 5;
wlen = 150;
zsc = 1.96; % 2.58 for 99% or 1.96 for 95%

% dummy_spks = sort(datasample(sess.ind(root.tsb(1):root.tsb(end)),15000,'Replace',false));
% cc = dummy_spks';

[orig_count,binedges] = plot_frXripple(root,cc,sess,wlen,histoBnsz,1);
tmpRipMod = get_RipMod(root,cc,sess,1000);
tmpRipPart = get_RipParticipation(root,cc,sess,wlen);
legend(['P(Participation) ' num2str(tmpRipPart)])
jit_count = zeros(nShufs,length(orig_count));

for i = 1:nShufs
    % Circular shuffle method
    shiftRoot = shiftTrain(root,sess);
    jit_count(i,:) = plot_frXripple(shiftRoot,cc,sess,150,histoBnsz,0);

    % % Jitter shuffle method
    % jits = jitter_spks(root.tsb(root.cl == cc), jitlen);   
    % % jits = jitter_spks(cc, jitlen);   
    % jit_count(i,:) = plot_frXripple(root,jits,sess,wlen,histoBnsz,0);
end

% for i = 1:length(root.good)
%     jit_count(i,:) = plot_frXripple(root,root.good(i),150,histoBnsz,0);
% end

% Try Fujisawa et al. 2008 point-wise and global confidence bands
% For all shuffles, find how often orig_count fell above or below
p_over = sum(jit_count >= orig_count) / (nShufs + 1);
p_undr = sum(jit_count <= orig_count) / (nShufs + 1);
sig_over = p_over < 0.025;   % alpha 0.05/2
sig_undr = p_undr < 0.025;

% Global confidence band is proportion of shuffled data sets that had even
% 1 location above/below the associated pointwise band
% Iteratively search stricter pointwise band alpha until global, Family-wise alpha = 5%
gbl_alpha = 1.0;
pt_delta = 0;
while gbl_alpha > 0.05
    tmpPctUp = prctile(jit_count,97.5+pt_delta,1);
    tmpPctDn = prctile(jit_count,2.5-pt_delta,1);
    gbl_cross = jit_count > tmpPctUp | jit_count < tmpPctDn;
    gbl_alpha = sum(sum(gbl_cross,2) > 0)/nShufs;
    pt_delta = pt_delta + 0.05;
    % Error out if get above 100 percentile
    if pt_delta + 97.5 > 100
        gbl_alpha = 0;
        pt_delta = pt_delta + 0.05; % Ensure final check gets to 100th pctile
    end
end
pt_delta = pt_delta - 0.05; %Correct for final loop iteration adjustment

% Final check if true data exceeds globally set pointwise bands at any point
p_final = orig_count > prctile(jit_count,97.5+pt_delta,1) | orig_count < prctile(jit_count,2.5-pt_delta,1);

figure; hold on
xcoords = binedges(1:end-1)+histoBnsz/2000;
% errorbar(xcoords,binMu,zsc*(binStd / sqrt(nShufs)),'k')
bar(xcoords,orig_count,'FaceColor','flat','CData',[.5 .5 .5]);
plot(xcoords,prctile(jit_count,97.5,1),'b-')
plot(xcoords,prctile(jit_count,2.5,1),'b-','HandleVisibility','off')
plot(xcoords,prctile(jit_count,97.5+pt_delta,1),'m-')
plot(xcoords,prctile(jit_count,2.5-pt_delta,1),'m-','HandleVisibility','off')
if sum(p_final) > 0
    plot(xcoords(sig_over | sig_undr),orig_count(sig_over | sig_undr)+1,'k*')
end
xlabel('Time to ripple peak (sec)'); ylabel('Spike Count')
legend({'Real','95% Pointwise','95% Global','Pt-wise < 0.05'})
set(gca,'FontSize',12,'FontName','Arial')

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
        [orig_count,binedges] = plot_frXripple(root,cc,sess,wlen,histoBnsz,0);
        ripMod(ct,1:3) = [cc, j-1 root.info.uType(lyrCCs(i))];
        ripMod(ct,4) = get_RipParticipation(root,cc,sess,wlen);

        jit_count = zeros(nShufs,length(orig_count));
        for k = 1:nShufs
            % Circular shuffle method
            shiftRoot = shiftTrain(root,sess);
            jit_count(k,:) = plot_frXripple(shiftRoot,cc,sess,wlen,histoBnsz,0);

            % % Jitter shuffle method
            % jits = jitter_spks(root.tsb(root.cl == cc), jitlen);
            % % jits = jitter_spks(cc, jitlen);
            % jit_count(i,:) = plot_frXripple(root,jits,sess,wlen,histoBnsz,0);
        end

        gbl_alpha = 1.0;
        pt_delta = 0;
        while gbl_alpha > 0.05
            tmpPctUp = prctile(jit_count,97.5+pt_delta,1);
            tmpPctDn = prctile(jit_count,2.5-pt_delta,1);
            gbl_cross = jit_count > tmpPctUp | jit_count < tmpPctDn;
            gbl_alpha = sum(sum(gbl_cross,2) > 0)/nShufs;
            pt_delta = pt_delta + 0.05;
            % Error out if get above 100 percentile
            if pt_delta + 97.5 > 100
                gbl_alpha = 0;
            end
        end
        pt_delta = pt_delta - 0.05; %Correct for final loop iteration adjustment

        % Final check if true data exceeds globally set pointwise bands at any point
        p_final = orig_count > prctile(jit_count,97.5+pt_delta,1) | orig_count < prctile(jit_count,2.5-pt_delta,1);

        ripMod(ct,5) = sum(p_final);    % If more than 1 point falls outside global significance
        ct = ct + 1;
    end
    disp(['Finished shuffling shank ' num2str(j)])
end

root.ripMod = ripMod;

%% Observational statistics

disp(['Overall in-layer significant ' num2str(sum(ripMod(:,5))) '/' num2str(size(ripMod,1))])

for i = 1:nshanks
    disp(['Shank ' num2str(i-1) ' significant ' num2str(sum(ripMod(ripMod(:,2) == i-1,5))) '/' num2str(numel(ripMod(ripMod(:,2) == i-1,5)))])
    disp(['      ' num2str(i-1) ' mean participation: ' num2str(mean(ripMod(ripMod(:,2) == i-1,4)))])
end

%% Waterfall
histoBnsz = 1;
dummy_counts = plot_frXripple(root,1,sess,wlen,histoBnsz,0);
gaussSigma = 5;
sz = 32;
gaussX = linspace(-sz/2, sz/2, sz);
gaussFilt = exp(-gaussX .^2 / (2*gaussSigma ^2));
gaussFilt = gaussFilt / sum(gaussFilt); % Normalize

ct = 1;
ripZMap = [];

for i = 1:size(ripMod,1)
    if ripMod(i,5) == 1
        cc = ripMod(i,1);
        [tmpPeriSpikes,bins] = plot_frXripple(root,cc,sess,wlen,histoBnsz,0);
        muRipSpks = mean(tmpPeriSpikes);
        sdRipSpks = std(tmpPeriSpikes);
        zRipSpks = (tmpPeriSpikes - muRipSpks) ./ sdRipSpks;
        ripZMap(ct,:) = conv(zRipSpks,gaussFilt,'same');
    else
        ripZMap(ct,:) = zeros(1,length(dummy_counts));
    end
    ct = ct + 1;
end

sigMod = ripMod(:,5) == 1;
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










