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
sess.name = bhvrFile.name;

% Load data from .h5 file
dat             = ws.loadDataFile(sess.name);
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
[~,sess.lckind] = findpeaks(double(sess.lck > 0.5));
sess.slx        = double(sess.aidat(:,8) > 0.5);    % Translate to binary

% Digital data; multiplexed from all digital channels; 1 = LapReset; 2 = Reward release
sess.didat      = dat.(datFields{2}).digitalScans;
sess.rwd        = double(sess.didat == 2);
[~,sess.rwdind] = findpeaks(sess.rwd);
sess.rst        = double(sess.didat == 1);
[~,sess.lapstt] = findpeaks(sess.rst); sess.lapstt = sess.lapstt + 1; % Account for position reset lagged by 1 index
sess.lapstt     = [sess.ind(1); sess.lapstt];   %use first ts as first lap start
sess.lapend     = [sess.lapstt(2:end) - 1; sess.ind(end)];     %Use last ts as last lap end
sess.nlaps      = size(sess.lapend,1);

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

end