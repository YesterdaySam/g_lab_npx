function [fhandle] = fixRatio(fhandle, xratio, yratio)
%% Converts aspect ratio from landscape to portrait
%
% Inputs:
%   fhandle = figure handle
%   xratio = 0.5625; divisor for x dimension of fhandle
%   yratio = 1.7778; divisor for y dimension of fhandle
%
% Outputs:
%   fhandle = updated figure (or original, if screen is still portrait)
%
% Created 8/3/2026 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments 
    fhandle
    xratio = 0.5625  % Convert 1920 to 1080
    yratio = 1.7778  % Convert 1080 to 1920
end

aspectRatio = get(0,'ScreenSize');

if aspectRatio(3) < aspectRatio(4)  % No change unless screen is in portrait mode
    newX = fhandle.Position(3) / xratio;
    newY = fhandle.Position(4) / yratio;
% else
    % newX = fhandle.Position(3) * xratio;
    % newY = fhandle.Position(4) * yratio;
end
set(gcf,'units','normalized','position',[fhandle.Position(1) fhandle.Position(2) newX newY])

end