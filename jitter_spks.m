function [jit_spks] = jitter_spks(orig_spks,wnlen)
%% Uniformly, independently jitter spikes on the interval [-wnlen wnlen]
%
% Inputs:
% orig_spks = spiketrain indices
% wnlen = interval in indices over which to jitter
%
% Outputs:
% jit_spks = jittered spike train indices
%
% Created 5/12/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    orig_spks
    wnlen = 25; % Assuming 2500Hz behavior sampling rate, +/- 10 msec
end

jitoffset = randi([-wnlen,wnlen],length(orig_spks),1);
jit_spks = orig_spks + jitoffset;

end