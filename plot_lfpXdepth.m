function [fhandle] = plot_lfpXdepth(root)
%% Plots estimated power spectral density by shank, depth and freq band
%
% Inputs:
% root = root object. Must have lyrbounds and uPSD fields
%
% Outputs:
% fhandle = handle to figure
%
% Created 2/15/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------


try
    disp(['Existing root.lyrbounds: ' num2str(root.lyrbounds(1,:)) ' TO ' num2str(root.lyrbounds(2,:))])
catch
    root = get_layerUnits(root,100);
    disp(['Added root.lyrbounds: ' num2str(root.lyrbounds(1,:)) ' TO ' num2str(root.lyrbounds(2,:))])
end

cmap(1).map = cool(4);
cmap(2).map = hot(8);
cmap(3).map = gray(5);

fhandle = figure; hold on
legcell = {};
nChXSh = height(root.lfpinfo)/numel(unique(root.lfpinfo.lfpShank));
grid on
for sh = 0:3
    for band = 1:size(root.bands,1)
        plot(root.uPSD(band,root.lfpinfo.lfpShank == sh),root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh),'Color',cmap(band).map(sh+1,:))
        legcell(band) = {[num2str(root.bands(band,1)) '-' num2str(root.bands(band,2)) 'Hz']};
    end
    for band = 1:size(root.bands,1)
        [tmpMax, tmpInd] = max(root.uPSD(band,root.lfpinfo.lfpShank == sh));
        tmpD = root.lfpinfo.lfpDepth(root.lfpinfo.lfpShank == sh);
        plot(tmpMax+1,tmpD(tmpInd),'<','Color',cmap(band).map(sh+1,:))
        if band == 2; plot([tmpMax tmpMax]+1,root.lyrbounds(:,sh+1),':_','Color',cmap(band).map(sh+1,:)); end
    end
end
xlabel('PSD (dB/Hz)')
ylabel('Distance from Probe tip (um)')
legend(legcell)
set(gca,'FontSize',12,'FontName','Arial')

end