function [dErr] = get_circDecodeErr(rPos,dPos,tracklen)
%% Get circularly computed decoding error between rPos and dPos
% If rPos = 0.1 and dPos = 1.8, and tracklen = 1.85, dErr = 0.06
%
% Inputs:
%   rPos     = Nx1 vector of real positions
%   dPos     = Nx1 vector of decoded positions
%   tracklen = double, in meters e.g. 1.85
%
% Outputs:
%   dErr = circular error between rPos and dPos
%
% Created 10/22/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    rPos
    dPos
    tracklen
end

wrap = tracklen/2;
rawErr = rPos - dPos;
dErr = rawErr;
dErr(rawErr > wrap) = dErr(rawErr > wrap) - tracklen;
dErr(rawErr < -wrap) = dErr(rawErr < -wrap) + tracklen;

end