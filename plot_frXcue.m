function [fhandle,rastfig] = plot_frXcue(root,unit,sess,qstruc,dbnsz,plotflag)
%% Plots the avg binned firing rate centered around times in qstruc
% For speed, pre-calculate random cue starts with get_QTwin.m function
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   unit = cluster ID
%   sess = session struct from importBhvr
%   qstruc = qs.q, qs.stt, qs.end = time indices; qs.qID = trial IDs
%   dbnsz = size of position bins, default 0.05m = 5cm
%   plotflag = binary of whether to plot the output
%
% Outputs:
%   binedges = spatial bin edges
%   binfr = spatial-binned firing rate
%   fhandle = handle to figure
%
% Created 1/27/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    unit {double}   %Cluster ID
    sess            %session struct
    qstruc          %struct with qs.q, qs.stt, qs.end, qs.qID
    dbnsz = 0.05    %m
    plotflag = 1    %binary
end

[~,~,~,~,spkmap,bnoccs,rastMt,binfr,binedges] = get_SI_cue(root,unit,sess,qstruc,dbnsz);

spkct = sum(spkmap,1);
occct = sum(bnoccs,1);

if plotflag
    rawfr = spkct ./ occct;
    rawfr(isnan(rawfr)) = 0;
    sem = std(spkmap ./ bnoccs,'omitnan')/sqrt(sess.nlaps);
    sem(isnan(sem)) = 0;
    ciup = rawfr + sem*1.96; %weird for some laps big SEM?
    cidn = rawfr - sem*1.96;
    
    fhandle = figure; hold on
    plot([binedges(1:length(rawfr))]*100,rawfr, 'k')
    patch(100*[binedges(1:length(rawfr)),fliplr(binedges(1:length(cidn)))],[cidn,fliplr(ciup)],'k','FaceAlpha',0.5,'EdgeColor','none')
    plot([binedges(1:end-1)]*100,binfr, 'b')

    xlabel('Distance from Cue (cm)'); 
    ylabel('Firing Rate (spk/s)')
    if max(mean(binfr,1,'omitnan'),[],'all') < 10
        ylim([0 10])
    elseif max(binfr,[],'all') < 20
        ylim([0 20])
    elseif max(binfr,[],'all') < 30
        ylim([0 30])
    end
    title(['Unit ' num2str(unit)])
    set(gca,'FontSize',12,'FontName','Arial')
    
    % Raster
    rastfig = figure; plot(rastMt(:,1)*100,rastMt(:,2),'k|')
    xlabel('Distance from Cue (cm)'); 
    ylabel('Trial')
    set(gca,'FontSize',12,'FontName','Arial')
end

end