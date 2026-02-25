function [root] = get_layerBounds(root,minWidth,minProm)
%% Assigns units based on peak in ripple band power
%
% Inputs:
% root = root object. Must have uPSDMax field
% minWidth = minimum layer width (100um) or use string to use peak-defined
% fields for all
%
% Outputs:
% root = updated root with root.lyrbounds
%
% Created 3/11/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    minWidth = 100  % In um
    minProm  = 0.25 % Percent of prominence at which to assign lyrBound
end
nShanks = numel(unique(root.lfpinfo.lfpShank));

for sh = 0:nShanks-1
    tmpsh_psd = root.uPSD(2,root.lfpinfo.lfpShank == sh);
    [tmpmax, tmpmaxind] = max(tmpsh_psd);
    tmpsh_depth = root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh)';

    if tmpmaxind == 1 || tmpmaxind == 2 %Pad on either side to enable find peaks on the edges
        tmpsh_psd = [mean(tmpsh_psd), tmpsh_psd]; 
        tmpsh_depth = [NaN, tmpsh_depth];
    elseif tmpmaxind == length(tmpsh_psd) || tmpmaxind == length(tmpsh_psd)-1
        tmpsh_psd = [tmpsh_psd, mean(tmpsh_psd)];
        tmpsh_depth = [tmpsh_depth, NaN];
    end

    [pks,locs,widths,proms] = findpeaks(tmpsh_psd,'MinPeakProminence',std(tmpsh_psd)); %Get peaks and prominences
    if isempty(pks) % In case of small double peak near mean peak leading to spurious prominence
        % [pks,locs,widths,proms] = findpeaks(smooth(tmpsh_psd,3),'MinPeakProminence',std(tmpsh_psd)*.5);
        [pks,locs,widths,proms] = findpeaks(smooth(tmpsh_psd,3));
    end
    tmpind = find(pks == max(pks));    %Index of highest peak
    tmploc = locs(tmpind); %Channel index of highest peak
    tmpprom = proms(tmpind);

    % Split PSD into two vectors and search for last point of minProm prominence in either direction to get layer boundaries
    v1 = fliplr(tmpsh_psd(1:tmploc));
    v2 = tmpsh_psd(tmploc:end);
    lim1 = find(diff(v1 > pks(tmpind) - tmpprom*(1-minProm)) == -1,1,'first');
    lim2 = find(diff(v2 > pks(tmpind) - tmpprom*(1-minProm)) == -1,1,'first');
    lowD = tmpsh_depth(tmpsh_psd == v1(lim1));
    uppD = tmpsh_depth(tmpsh_psd == v2(lim2));
    
    root.lyrbounds(:,sh+1) = [lowD; uppD];
end

%% Handle narrow layer bounds by padding out up to half minWidth, but don't go beyond 0 or max depth of probe

if isnumeric(minWidth)
    failMinWidth = (root.lyrbounds(2,:) - root.lyrbounds(1,:)) < minWidth;
    for i = 1:nShanks
        if failMinWidth(i)
            extLow = root.lfpinfo.lfpDepth(root.uPSDMax(2,i)) - minWidth/2;
            if extLow < root.lyrbounds(1,i) && extLow >= 0
                root.lyrbounds(1,i) = extLow;
            elseif extLow < root.lyrbounds(1,i) %Edge case of low bound below 0
                root.lyrbounds(1,i) = 0;
            end
            extHi = root.lfpinfo.lfpDepth(root.uPSDMax(2,i)) + minWidth/2;
            if extHi > root.lyrbounds(2,i) && extHi <= max(root.info.depth)
                root.lyrbounds(2,i) = extHi;
            elseif extHi > root.lyrbounds(2,i)
                root.lyrbounds(2,i) = max(root.info.depth);
            end
        end
    end
end

end