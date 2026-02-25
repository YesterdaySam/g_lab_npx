function [decodeInfo] = decodePosBayes(root,sess,tau,useUnits)
%%% Use

arguments
    root
    sess
    tau = 0.5   % Time Window, seconds
    useUnits = root.goodind     % Binary vector 1xN where N = all units
end

decodeInfo = [];
newts = sess.ts(sess.runInds & sess.lapInclude);    % Use only run periods
newpos = sess.pos(sess.runInds & sess.lapInclude);
dbnsz = 0.05;
if islogical(useUnits)
    ccs = root.info.cluster_id(useUnits);
else
    ccs = useUnits;
end

subsamp = sess.samprate / 50;
for i = 1:length(ccs)
    cc = ccs(i);
    [~,~,~,~,~,~,posfr(i,:),binedges] = get_SI(root,cc,sess,dbnsz);
end
binpos = binedges(1:end-1)+dbnsz/2;

expectSpk = posfr' + (eps.^8); % pos x cell
expectSpk = expectSpk * tau;
ct = 1;

for i = 1:subsamp:length(newts)/10
    if newts(i) - tau/2 <= 0 || newts(i) + tau/2 >= sess.ts(end)  % Ignore times before/after the minimum window
        continue
    end
    firstInd = find(sess.ts > newts(i)-tau/2,1);
    lastInd = find(sess.ts < newts(i)+tau/2,1,'last');
    realpos = mean(sess.pos(firstInd:lastInd));

    spks = root.ts > newts(i) - tau/2 & root.ts < newts(i) + tau/2;
    spkIds = root.cl(spks);
    nSpks = histcounts(spkIds,0:max(root.good)+1);  % Don't use groupcounts - about 2x slower!
    curSpk = nSpks(ccs);   % Spike counts in window for good units only
    % curSpk = curSpk(bothSIUnits);
    useTmp = curSpk>=0;     % Use all or use only those which spiked

    % Group counts attempt
    % tmpSpk = groupcounts(spkIds,0:max(root.good)+1,'IncludeEmptyGroups',true);
    % curSpk = tmpSpk(root.good)';   % Spike counts in window for good units only
    % curSpk = curSpk(bothSIUnits);
    % useTmp = curSpk > 0;

    % Bayes rule, decode current location
    tmp = bsxfun(@power, expectSpk(:,useTmp), curSpk(useTmp)); % [nPos x nTbin x nCell]
    tmp = prod(tmp,2);
    expon = exp(-sum(expectSpk(:,useTmp),2));     % Sum rate map for 1:N cells
    post = bsxfun(@times, tmp, expon);
    post = post./sum(post); % Normalization
    [~,id] = max(post); % decoded position is the one with max posterior prob
    decodeInfo(ct,1) = binpos(id);
    decodeInfo(ct,2) = realpos;
    decodeInfo(ct,3) = binpos(id) - realpos;
    ct = ct+1;
end

figure; plot(decodeInfo(:,2))
hold on; plot(decodeInfo(:,1))
