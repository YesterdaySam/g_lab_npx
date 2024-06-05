function [sess] = importBhvr(datdir)
%%
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
sess.vel        = sess.aidat(:,1) - min(sess.aidat(:,1));
sess.pos        = sess.aidat(:,2);
sess.lck        = sess.aidat(:,3);
[~,sess.lckind] = findpeaks(double(sess.lck > 0.5));
sess.slx        = sess.aidat(:,8);

% Digital data; multiplexed from all digital channels; 1 = LapReset; 2 = Reward release
sess.didat      = dat.(datFields{2}).digitalScans;
sess.rwd        = double(sess.didat == 2);
[~,sess.rwdind] = findpeaks(sess.rwd);
sess.rst        = double(sess.didat == 1);
[~,sess.lapstt] = findpeaks(sess.rst); sess.lapstt = sess.lapstt + 1; % Account for position reset lagged by 1 index
sess.lapend     = sess.lapstt - 1;
% sess.posend     = find(diff(sess.pos) < -0.3);
% sess.posstt     = sess.posend + 1;
sess.nlaps      = size(sess.lapend,1);

% %% Temp Visualization 
% tmpmax = max(session.aidat);
% tmpmax = repmat(tmpmax,[max(session.ind),1]);
% figure; hold on;
% for i = 1:8
%     plot(session.ts, i-1 + session.aidat(:,i)./tmpmax(:,i))
% end
