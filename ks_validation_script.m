% Practice Script for KS data import
addpath('D:\code')

ks4path = 'D:\Test\short_data_test\kilosort4'; % Default batchsize is 2.000
ks25path = 'D:\Test\KS25_test\kilosort25_2'; % Default batchsize 2.1824
root4 = loadKS(ks4path);
root2 = loadKS(ks25path);
recLen = 90;

%% Spike Hole visualization - Spike time derivative method
el = 125;       % Pick highest FR electrode (125 for KS4, 690 for KS25_2
batch_sz2 = 2.1824;
batch_sz4 = 2;

% el_samps = root.ts(root.cl == el);      % All spiks of that electrode
root2samps = root2.ts;                     % Try with all spikes
root2deriv = diff(root2samps);              % Time difference between spikes
root2lag = root2deriv > 0.007;             % Spike holes > 7msec
root2holes = diff(root2samps(root2lag));    % Time difference between spikes on either side of spike hole
root2norm = root2holes/batch_sz2;           % Normalize to batch duration of KS

root4norm = diff(root4.ts(diff(root4.ts) > 0.007))/batch_sz4;

figure; histogram(root2norm,0:0.02:10)
hold on; histogram(root4norm,0:0.02:10)
ylabel('Spike Holes > 7ms (count)')
xlabel('Batch')
legend({'KS2.5','KS4'})

%% Spike Hole visualization - Binned FR method across whole probe

bnedges = 0:0.005:recLen; % ms bins
spkcts4 = histcounts(root4.ts,bnedges);
spkcts2 = histcounts(root2.ts,bnedges);

figure; plot(bnedges(1:end-1),spkcts2)
hold on; plot(bnedges(1:end-1),spkcts4)
plot([0 recLen],[0 0],'k--')
ylim([-5 inf])
ylabel('Total spikes (binned at 5ms')
xlabel('Time (s)')
legend({'KS2.5','KS4'})

%%  Spike Hole visualization - batch-end-aligned edges
msec = 0.005;
bnedges2 = 0:batch_sz2*msec:recLen;
bnedges4 = 0:batch_sz4*msec:recLen;

spkcts4 = histcounts(root4.ts,bnedges4);
spkcts2 = histcounts(root2.ts,bnedges2);

spk2ln = size(spkcts2,2);
spk2rm = mod(spk2ln,1/msec);

spkbatch2 = reshape(spkcts2(1:spk2ln-spk2rm), [1/msec,(spk2ln-spk2rm)*msec]);
spkbatch4 = reshape(spkcts4, [1/msec, size(spkcts4,2)*msec]);

figure; plot(linspace(0,batch_sz2,1/msec),sum(spkbatch2'))
hold on; plot(linspace(0,batch_sz4,1/msec),sum(spkbatch4'))
ylabel('Total batched spikes (binned at 5ms')
xlabel('Time (msec)')
legend({'KS2.5','KS4'})

%% Prep roots with basic metrics
nc2 = size(unique(root2.cl),1);
nc4 = size(unique(root4.cl),1);
root2.id = unique(root2.cl);
root4.id = unique(root4.cl);

for i = 1:nc2
    root2.nspks(i) = size(root2.ts(root2.cl == root2.id(i)),1);
end
root2.fr = root2.nspks/recLen;

for i = 1:nc4
    root4.nspks(i) = size(root4.ts(root4.cl == root4.id(i)),1);
end
root4.fr = root4.nspks/recLen;

%% Compare well-matched units with different metrics
matchT = readtable('D:\Test\short_data_test\KS_compare_units.xlsx');

ks4units = matchT.KS4;
ks2units = matchT.KS2;

for i = 1:length(ks4units)
    mapper4(i) = find(root4.id == ks4units(i));
    mapper2(i) = find(root2.id == ks2units(i));
end

ps = struct;
stats = struct;
corrType = 'Spearman';


[stats.fr_rho, ps.fr_corr] = corr(root4.fr(mapper4)',root2.fr(mapper2)',"Type",corrType);
[stats.nspk_rho, ps.nspk_corr] = corr(root4.nspks(mapper4)',root2.nspks(mapper2)',"Type",corrType);

figure; hold on
scatter(root2.fr(mapper2),root4.fr(mapper4))
plot([0 max(root4.fr)], [0 max(root4.fr)], 'k--')
xlabel('KS2.5 FR')
ylabel('KS4 FR')

figure; hold on
scatter(root2.nspks(mapper2),root4.nspks(mapper4))
plot([0 max(root4.nspks)], [0 max(root4.nspks)], 'k--')
xlabel('KS2.5 # Spikes')
ylabel('KS4 # Spikes')

%%
[~,ps.fr_tt,~,stats.fr_tt] = ttest(root2.nspks(mapper2) - root4.nspks(mapper4));

figure; hold on
histogram(root2.nspks(mapper2) - root4.nspks(mapper4),50)
plot([0 0], [0 20], 'k--','LineWidth',1)
xlabel('\Delta # Spikes KS2.5 - KS4')
ylabel('Count')

%% Test SI KS4 vs Raw KS4
datpath2 = 'D:\Kelton\probe_data\KW005\KW005_07182024_rec_D3_CA1_g0\KW005_07182024_rec_D3_CA1_g0_imec0';
datpath1 = 'D:\Kelton\probe_data\KW005\KW005_07182024_rec_D3_CA1_g0\preprocess';
kssipath = 'D:\Kelton\probe_data\KW005\KW005_07182024_rec_D3_CA1_g0\preprocess\kilosort4_direct'; % Default batchsize is 2.000
ksrawpath = 'D:\Kelton\probe_data\KW005\KW005_07182024_rec_D3_CA1_g0\KW005_07182024_rec_D3_CA1_g0_imec0\kilosort4'; % Default batchsize 2.1824
root_si = loadKS(datpath1,kssipath);
root_raw = loadKS(datpath2,ksrawpath);
recLen = root_raw.tspulse(end);

%% Add basic metrics to each unit
nc_raw = size(unique(root_raw.cl),1);
nc_si = size(unique(root_si.cl),1);
root_raw.id = unique(root_raw.cl);
root_si.id = unique(root_si.cl);

for i = 1:nc_raw
    root_raw.nspks(i) = size(root_raw.ts(root_raw.cl == root_raw.id(i)),1);
end
root_raw.fr = root_raw.nspks/recLen;

for i = 1:nc_si
    root_si.nspks(i) = size(root_si.ts(root_si.cl == root_si.id(i)),1);
end
root_si.fr = root_si.nspks/recLen;

%% Compare well-matched units with different metrics
matchT = readtable('D:\Kelton\probe_data\KW005\KW005_07182024_rec_D3_CA1_g0\preprocess\KS_compare_units.xlsx');
matchT2 = rmmissing(matchT(:,1:4)); %Remove unmatched units

rawunits = matchT2.KS4_raw;
siunits = matchT2.KS4_SI_direct;

for i = 1:length(rawunits)
    mapper_raw(i) = find(root_raw.id == rawunits(i));
    mapper_si(i) = find(root_si.id == siunits(i));
end

ps = struct;
stats = struct;
corrType = 'Spearman';


[stats.fr_rho, ps.fr_corr] = corr(root_raw.fr(mapper_raw)',root_si.fr(mapper_si)',"Type",corrType);
[stats.nspk_rho, ps.nspk_corr] = corr(root_raw.nspks(mapper_raw)',root_si.nspks(mapper_si)',"Type",corrType);

figure; hold on
scatter(root_si.fr(mapper_si),root_raw.fr(mapper_raw))
plot([0 max(root_raw.fr)], [0 max(root_raw.fr)], 'k--')
xlabel('KS SI FR')
ylabel('KS Raw FR')

figure; hold on
scatter(root_si.nspks(mapper_si),root_raw.nspks(mapper_raw))
plot([0 max(root_raw.nspks)], [0 max(root_raw.nspks)], 'k--')
xlabel('KS SI # Spikes')
ylabel('KS Raw # Spikes')

%%
[~,ps.fr_tt,~,stats.fr_tt] = ttest(root_si.nspks(mapper_si) - root_raw.nspks(mapper_raw));

figure; hold on
histogram(root_si.nspks(mapper_si) - root_raw.nspks(mapper_raw),50)
plot([0 0], [0 50], 'k--','LineWidth',1)
xlabel('\Delta # Spikes KS SI - Raw')
ylabel('Count')
