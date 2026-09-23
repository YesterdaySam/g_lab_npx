function [decodeI] = decodePosBayes(root,sess,expectSpk,useUnits,tau)
%% Use Bayes' rule to decode position over time from useUnits
% Subsamples behavior to 50Hz and estimates position from spikes within tau
% using the prior probability estimate from expectSpk (mean FR mat)
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   sess = session struct from importBhvr
%   expectSpk = NxM matrix of firing rates of N units by M positions
%   useUnits    % Inds of units e.g. root.good
%   tau = 0.5   % Time Window within which to count spikes, seconds
%
% Outputs:
%   decodeI = struct of decoded information
%       rPos = real absolute track position (m)
%       dPos = estiamted absolute track position (m)
%       dMat = NxT matrix of N posterior estimates at T time points
%       dErr = Circularly corrected error between rPos and dPos (m)
%       newT = 50Hz subsampled time stamps 
%
% Created 9/8/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    sess
    expectSpk   % NxM = Firing rate of N units by M positions
    useUnits    % Inds of units e.g. root.good
    tau = 0.5   % Time Window, seconds
end

newts = sess.ts(sess.runInds & sess.lapInclude);    % Use only run periods
newpos = sess.pos(sess.runInds & sess.lapInclude);
dbnsz = 0.05;

binpos = dbnsz/2:dbnsz:sess.maxPos;

expectSpk = (expectSpk' + (eps.^8)) * tau; % M pos x N units

sumExpectSpk = sum(expectSpk, 2); % Sum of expected spikes across ALL included units for each position: M pos x 1
logExpectSpk = log(expectSpk);  % Pre calculate log expectation of M pos x N units

subsamp = sess.samprate / 50;
useInds = 1:subsamp:length(newts);
nTs = length(useInds);

dPos = zeros(nTs,1);
dMat = zeros(length(sumExpectSpk),nTs);

ct = 1;

for i = useInds
    t_cur = newts(i);
    t_stt = t_cur - tau/2;
    t_end = t_cur + tau/2;

    if t_stt <= 0 || t_end >= sess.ts(end)  % Ignore times before/after the minimum window
        continue
    end

    spks = root.ts > t_stt & root.ts < t_end;
    spkIds = root.cl(spks)+1;   % Account for 0-indexing
    nSpks = histcounts(spkIds,0.5:1:max(useUnits+1)+0.5); % Over all units, accounting for 0-indexing
    curSpk = nSpks(useUnits+1);   % only spikes from useUnits, accounting for 0-indexing

    % Bayesian decoding
    log_post = (logExpectSpk * curSpk') - sumExpectSpk;

    % Subtract max for numerical stability before exponentiating
    post = exp(log_post - max(log_post));
    post = post / sum(post); % Normalize to probability distribution

    % % Bayes rule, decode current location
    % tmp = bsxfun(@power, expectSpk(:,useTmp), curSpk(useTmp)); % [nPos x nTbin x nCell]
    % tmp = prod(tmp,2);
    % expon = exp(-sum(expectSpk(:,useTmp),2));     % Sum rate map for 1:N cells
    % post = bsxfun(@times, tmp, expon);
    % post = post./sum(post); % Normalization

    [~,id] = max(post); % decoded position is the one with max posterior prob

    dPos(ct)    = binpos(id);
    dMat(:,ct)  = post;
    ct = ct+1;
end

rPos = newpos(useInds);

% Circularly calculate error
wrap = sess.maxPos/2;
rawErr = rPos - dPos;
dErr = rawErr;
dErr(rawErr > wrap) = dErr(rawErr > wrap) - sess.maxPos;
dErr(rawErr < -wrap) = dErr(rawErr < -wrap) + sess.maxPos;

% Remove unused preallocated rows and assign to output variable
validFrames  = 1:(ct-1);
decodeI.rPos = rPos(validFrames);
decodeI.dPos = dPos(validFrames);
decodeI.dMat = dMat(:,validFrames);
decodeI.dErr = dErr(validFrames);
decodeI.newT = newts(useInds(validFrames))';

end