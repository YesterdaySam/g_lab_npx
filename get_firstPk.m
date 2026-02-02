function [firstPkT] = get_firstPk(datStruc,unit)
%% Find first peak of FR relative to opto pulses for unit
%
% Inputs:
%   datStruc = datStruc object with optoMat and optobins fields
%   unit = cluster ID
%
% Outputs:
%   firstPkT = time (ms) of first post-pulse peak, or NaN if no peak
%
% Created 1/21/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------


tmpOpto = squeeze(datStruc.optoMat(unit,:,:));

baseInds = datStruc.optobins < 0;
counts = sum(tmpOpto);
baseMu = mean(counts(baseInds));
baseSD = std(counts(baseInds));

countZ = (counts - baseMu) ./ baseSD;
smoothZ = smoothdata(countZ,'gaussian',10);

[tmppks, tmplocs] = findpeaks(smoothZ);

% figure; hold on
% plot(datStruc.optobins,countZ)
% plot(datStruc.optobins,smoothZ)
% plot(datStruc.optobins(tmplocs),tmppks,'v')

% first peak >2 Z score
postPulsePks = datStruc.optobins(tmplocs) > 0;
hiPks = tmppks > 2;

firstPkT = datStruc.optobins(tmplocs(find(postPulsePks & hiPks, 1)));

if isempty(firstPkT)
    firstPkT = NaN;
end

end
