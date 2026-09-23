function [pvCorrXTime] = get_pvXtime(refMap,frMap,useUnits,useMask)
%% Returns the population vector of frMap against refMap
%
% Inputs:
% refMap = MxN firing rate of M units over N spatial bins 
% frMap = LxNxM firing rate of M units over N bins over L laps
% useUnits = binary of Mx1 of units to include
% useMask = binary of NxN positions to use, otherwise use whole corr
%
% Outputs:
% pvCorrXTime = LxN Pop. Vector correlation in NxN space over L laps
%
% Created 8/6/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    refMap
    frMap
    useUnits = 0
    useMask = 0
end

if useUnits == 0
    useUnits = ones(size(frMap,3));
end

if useMask == 0
    useMask = ones(size(frMap,2));
end

% Ensure refMap is normalized [0 1]
refMap = refMap(useUnits,:);
refMap = refMap ./ max(refMap,[],2);

pvCorrXTime = zeros(size(frMap,1),1);

for i = 1:size(frMap,1)
    lapMap = squeeze(frMap(i,:,useUnits))';
    normFR = lapMap ./ max(lapMap,[],2);
    % tmpFR = normFR;
    % tmpFR(isnan(normFR(:,1)),:) = [];
    % plot_unitWaterfall(tmpFR,0:0.05:1.85);
    corrTmp = corr(normFR,refMap,'rows','complete');
    % plot_pvcorr(corrTmp.*useMask);
    % plot_pvcorr(corrTmp);
    pvCorrXTime(i) = mean(corrTmp(useMask),'all');
end

end