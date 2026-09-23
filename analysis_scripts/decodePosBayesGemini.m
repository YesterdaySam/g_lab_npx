function [decodeInfo] = decodePosBayesGemini(root, sess, useUnits, tau)
arguments
    root
    sess
    useUnits = root.goodind     % Binary vector 1xN or array of cluster IDs
    tau (1,1) double = 0.5      % Time Window in seconds
end

% Extract time and position during valid running laps
newts = sess.ts(sess.runInds & sess.lapInclude);
dbnsz = 0.05; % Position bin size

if islogical(useUnits)
    ccs = root.info.cluster_id(useUnits);
else
    ccs = useUnits;
end

% 1. Compute spatial tuning curves / rate maps (posfr: nUnits x nPos)
nUnits = length(ccs);
for i = 1:nUnits
    cc = ccs(i);
    [~,~,~,~,~,~, posfr(i,:), binedges] = get_SI(root, cc, sess, dbnsz);
end

binpos = binedges(1:end-1) + dbnsz/2;
nPosBins = length(binpos);

% Rate map matrix: [nPos x nUnits]
% Add epsilon to prevent log(0)
expectSpk = (posfr' + 1e-10) * tau; 

% Sum of expected spikes across ALL included units for each position: [nPos x 1]
sumExpectSpk = sum(expectSpk, 2); 
logExpectSpk = log(expectSpk);

subsamp = round(sess.samprate / 50);
validIndices = 1:subsamp:length(newts);

% Preallocate output matrix
decodeInfo = nan(length(validIndices), 3);
ct = 1;

% Pre-filter spike data range
maxClusterID = max(ccs);

for i = validIndices
    t_center = newts(i);
    t_start  = t_center - tau/2;
    t_end    = t_center + tau/2;
    
    if t_start <= 0 || t_end >= sess.ts(end)
        continue
    end
    
    % Interpolate or extract ground truth position
    firstInd = find(sess.ts > t_start, 1);
    lastInd  = find(sess.ts < t_end, 1, 'last');
    if isempty(firstInd) || isempty(lastInd) || firstInd > lastInd
        continue;
    end
    realpos = mean(sess.pos(firstInd:lastInd));
    
    % Extract spikes in current window
    spkIdx = root.ts > t_start & root.ts < t_end;
    spkIds = root.cl(spkIdx);
    
    % Bin spikes cleanly matching unit cluster IDs:
    % Bin 1 catches ID=1 ([0.5, 1.5]), Bin k catches ID=k ([k-0.5, k+0.5])
    counts = histcounts(spkIds, 0.5 : 1 : (maxClusterID + 0.5));
    curSpk = counts(ccs); % Spike count vector for specified units: [1 x nUnits]
    
    % 2. Log-Bayesian Decoding
    % log P(x|N) ~ sum(n_i * log(lambda_i(x))) - sum(lambda_i(x))
    log_post = (logExpectSpk * curSpk') - sumExpectSpk;
    
    % Subtract max for numerical stability before exponentiating
    post = exp(log_post - max(log_post)); 
    post = post / sum(post); % Normalize to probability distribution
    
    % Decoded position (MAP estimate)
    [~, id] = max(post);
    
    decodeInfo(ct, 1) = binpos(id);
    decodeInfo(ct, 2) = realpos;
    decodeInfo(ct, 3) = binpos(id) - realpos;
    ct = ct + 1;
end

% Remove unused preallocated rows
decodeInfo(ct:end, :) = [];

% Visualization
figure;
plot(decodeInfo(:, 2), 'k-', 'DisplayName', 'Actual Position');
hold on;
plot(decodeInfo(:, 1), 'r.', 'DisplayName', 'Decoded Position');
xlabel('Time Bins');
ylabel('Position');
legend();
title('Bayesian Position Decoding');