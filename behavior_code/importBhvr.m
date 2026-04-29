function [sess] = importBhvr(datdir)
%% Create a session variable for organizing neural and behavioral data
% Requires ws.loadDataFile functions from https://wavesurfer.janelia.org/
%%% Inputs
% datdir = filepath to directory containing an '.h5' file from wavesurfer
%%% Outputs
% sess = struct containing the following information
%   header      = wavesurfer header information
%   samprate    = behavioral sampling rate (e.g. 2KHz)
%   aidat       = Mx8 matrix of analog channel data, where M is all samples (see comments below)
%   ind         = indices of all
%   ts          = time stamps, spanning time of session
%   vel         = speed information, in m/s (may be translated by fixed value
%   pos         = position information, in m
%   lck         = lick port state from ~0 to ~10
%   lckind      = indices of lick onsets
%   slx         = up and down states from Imec card at 1Hz, 500ms Pulse Width
%   didat       = Mx1 array of multiplexed digital channel data, where M is all samples (see comments below)
%   rwd         = Mx1 array of just the digital reset signal ('2's)
%   rwdind      = indices of reward pulse onsets
%   rst         = Mx1 array of just the digital reset signal ('1's)
%   lapstt      = indices of lap starts
%   lapend      = indices of lap ends
%   nlaps       = number of laps (based on number of reset pulses)
%
% Created 5/15/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

clear sess dat datFields
% Identify .h5 data file
sess.bhvrDir = datdir;
cd(sess.bhvrDir)
bhvrFile = dir("*.h5");
sess.name = bhvrFile(1).name;

for i = 1:size(bhvrFile,1)
    % Load data from .h5 file
    dat             = ws.loadDataFile(bhvrFile(i).name);
    datFields       = fieldnames(dat);      % Use in case sweep number is not 001
    sess.header     = dat.header;
    sess.samprate   = sess.header.AcquisitionSampleRate;

    % Analog data: 1 = Vel; 2 = Dist (m); 3 = Licks; 4 = Reward; 5 = VmOut; 6 = Im; 7 = Bpod Target; 8 = SGLX
    sess.aidat      = dat.(datFields{2}).analogScans;
    sess.ind        = 1:size(sess.aidat,1);
    sess.ts         = sess.ind/sess.samprate;
    sess.vel        = sess.aidat(:,1);
    sess.velshft    = (sess.vel - mode(sess.vel))*80; % Subtract off velocity offset and convert analog voltage to cm/s, may leave negative values
    sess.pos        = sess.aidat(:,2);
    sess.lck        = sess.aidat(:,3);
    sess.opto       = sess.aidat(:,4);
    [~,sess.optoind]= findpeaks(double(sess.opto > 0.5));
    [~,sess.lckind] = findpeaks(double(sess.lck > 0.5));
    try
        sess.slx        = double(sess.aidat(:,8) > 0.5);    % Translate to binary
    catch
        warning('Failed to import synch pulse (sess.slx)')
    end

    % % for fox cue
    % try
        sess.randCue   = sess.aidat(:,5); %take care to account for switching
        sess.fixedCue    = sess.aidat(:,6);

        [~, sess.fixedCueInd] = findpeaks(double(sess.fixedCue > 0.5));
        [~, sess.randCueInd] = findpeaks(double(sess.randCue > 0.5));
    % 
    %     if ~isempty(sess.fixedCueInd) & ~isempty(sess.randCueInd) %if lights exist, within 500ms dont count
    %         keepFix = 1;
    %         keepRand = 1;
    %         for i = 2:length(sess.fixedCueInd)
    %             if sess.fixedCueInd(i) - sess.fixedCueInd(i-1) > 0.3*sess.samprate
    %                 keepFix(end+1) = i;
    %             end
    %         end
    % 
    %         for i = 2:length(sess.randCueInd)
    %             if sess.randCueInd(i) - sess.randCueInd(i-1) > 0.3*sess.samprate
    %                 keepRand(end+1) = i;
    %             end
    %         end
    %         sess.fixedCueInd = sess.fixedCueInd(keepFix);
    %         sess.randCueInd = sess.randCueInd(keepRand);
    %     end
    % catch
    % end

    % Digital data; multiplexed from all digital channels; 1 = LapReset; 2 = Reward release
    sess.didat      = dat.(datFields{2}).digitalScans;
    sess.rwd        = double(sess.didat == 2);
    [~,sess.rwdind] = findpeaks(sess.rwd);
    sess.rst        = double(sess.didat == 1);

    [~,sess.lapstt] = findpeaks(sess.rst);
    
    % Account for position holding steady for a few timestamps after reset
    for j = 1:length(sess.lapstt)
        tmp = find(sess.pos(sess.lapstt(j):end) < 0.1,1,'first');
        if tmp > 10; tmp = 1; disp('Extra long post-reset delay at track end in importBhvr'); end
        sess.lapstt(j) = sess.lapstt(j) + tmp -1;
    end
    % sess.lapstt     = sess.lapstt + 1;  % Account for position reset lagged by 1 index
    sess.lapstt     = [sess.ind(1); sess.lapstt];   %use first ts as first lap start
    sess.lapend     = [sess.lapstt(2:end) - 1; sess.ind(end)];     %Use last ts as last lap end

    sess.nlaps      = size(sess.lapend,1);
    sess.maxPos     = 1.85; % Hard coded to keep consistency
    % sess.maxPos     = round(median(sess.pos(sess.lapend)),2); % Track length based on lap ends
    sess            = getErrorTrials(sess);     % Identify trials of the right length and rewarded trials
    sess            = get_RunInds(sess,0.025,2); % Add runInds variable with binary of running (1) or standing (0)

    if sess.nlaps == 1      % In case of reset error
        try
            posrst = find(diff(sess.pos) < -0.3);
            sess.lapstt = [sess.ind(1); posrst+1];
            sess.lapend = [sess.lapstt(2:end) - 1; sess.ind(end)];
            sess.nlaps  = size(sess.lapend,1);
        catch
            warning(['Only 1 lap found for session ' sess.name])
        end
    end

    sess            = get_lapInclude(sess);

    sessOrg(i).sess = sess;
end

% Attempt to merge multiple sess objects
sess = sessOrg(1).sess;
if length(sessOrg) > 1
    for i = 2:length(sessOrg)
        iEnd = sess.ind(end);
        sess.ind        = [sess.ind, iEnd + (sessOrg(i).sess.ind)];
        sess.ts         = [sess.ts, sess.ts(end) + (sessOrg(i).sess.ts)];
        sess.velshft    = [sess.velshft; sessOrg(i).sess.velshft];
        sess.pos        = [sess.pos; sessOrg(i).sess.pos];
        sess.lck        = [sess.lck; sessOrg(i).sess.lck];
        sess.opto       = [sess.opto; sessOrg(i).sess.opto];
        sess.optoind    = [sess.optoind; iEnd + (sessOrg(i).sess.optoind)];
        sess.lckind     = [sess.lckind; iEnd + (sessOrg(i).sess.lckind)];
        sess.rwd        = [sess.rwd; sessOrg(i).sess.rwd];
        sess.rwdind     = [sess.rwdind; iEnd + (sessOrg(i).sess.rwdind)];
        sess.rst        = [sess.rst; sessOrg(i).sess.rst];
        sess.lapstt     = [sess.lapstt; iEnd + (sessOrg(i).sess.lapstt)];
        sess.lapend     = [sess.lapend; iEnd + (sessOrg(i).sess.lapend)];
        sess.errTrials  = [sess.errTrials, sess.nlaps + (sessOrg(i).sess.errTrials)];
        sess.valTrials  = [sess.valTrials, sess.nlaps + (sessOrg(i).sess.valTrials)];
        sess.rwdTrials  = [sess.rwdTrials, sess.nlaps + (sessOrg(i).sess.rwdTrials)];
        sess.nlaps      = sess.nlaps + sessOrg(i).sess.nlaps;
        sess.runInds    = [sess.runInds; sessOrg(i).sess.runInds];
        sess.lapInclude = [sess.lapInclude; sessOrg(i).sess.lapInclude];
    end
end

% % Use only the session with highest laps
% for i = 1:size(bhvrFile,1)
%     nlaps(i) = sessOrg(i).sess.nlaps;
% end
% [~, maxSess] = max(nlaps);
% sess = sessOrg(maxSess).sess;

% Remove extraneous fields
sess = rmfield(sess,'aidat');
sess = rmfield(sess,'didat');

end