function [p_final,fhandle] = get_confband(shufMat,orig_dat,plotFlag,binedges,bnsz)
%% Returns the items in orig_dat that fall beyond global confidence band
% Employs Fujisawa et al. 2008 point-wise and global confidence bands
% For all shuffles, find how often orig_count fell above or below
%
% Inputs:
% shufMat = MxN matrix of shuffled data, M = nShuffles, N = length of
%   orig_dat
% orig_dat = 1xN vector of original data values
%
% Outputs:
% p_final = spatial bin edges
% binfr = spatial-binned firing rate
% fhandle = handle to figure
%
% Created 7/15/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    shufMat
    orig_dat
    plotFlag = 0
    binedges = 1:length(orig_dat)
    bnsz = 0.5
end

nShufs = size(shufMat,1);

gbl_alpha = 1.0;
pt_delta = 0;

% Global confidence band is proportion of shuffled data sets that had even
% 1 location above/below the associated pointwise band
% Iteratively search stricter pointwise band alpha until global, Family-wise alpha = 5%

while gbl_alpha > 0.05
    tmpPctUp = prctile(shufMat,97.5+pt_delta,1);
    tmpPctDn = prctile(shufMat,2.5-pt_delta,1);
    gbl_cross = shufMat > tmpPctUp | shufMat < tmpPctDn;
    gbl_alpha = sum(sum(gbl_cross,2) > 0)/nShufs;
    pt_delta = pt_delta + 0.05;
    % Error out if get above 100 percentile
    if pt_delta + 97.5 > 100
        gbl_alpha = 0;
    end
end
pt_delta = pt_delta - 0.05; %Correct for final loop iteration adjustment

% Final check if true data exceeds globally set pointwise bands at any point
p_final = orig_dat > prctile(shufMat,97.5+pt_delta,1) | orig_dat < prctile(shufMat,2.5-pt_delta,1);

if plotFlag
    fhandle = figure; hold on
    xcoords = binedges+bnsz*0.5;
    p_over = sum(shufMat >= orig_dat) / (nShufs + 1);
    p_undr = sum(shufMat <= orig_dat) / (nShufs + 1);
    sig_over = p_over < 0.025;   % alpha 0.05/2
    sig_undr = p_undr < 0.025;

    bar(xcoords,orig_dat ./ sum(orig_dat,'all'),'FaceColor','flat','CData',[.5 .5 .5]);
    plot(xcoords,prctile(shufMat,97.5,1) ./ sum(orig_dat,'all'),'b-')
    plot(xcoords,prctile(shufMat,2.5,1) ./ sum(orig_dat,'all'),'b-','HandleVisibility','off')
    plot(xcoords,prctile(shufMat,97.5+pt_delta,1) ./ sum(orig_dat,'all'),'m-')
    plot(xcoords,prctile(shufMat,2.5-pt_delta,1) ./ sum(orig_dat,'all'),'m-','HandleVisibility','off')
    if sum(p_final) > 0
        plot(xcoords(sig_over | sig_undr),orig_dat(sig_over | sig_undr) ./ sum(orig_dat,'all')+.01,'k*')
    end
    ylabel('Probability')
    legend({'Real','95% Pointwise','95% Global','Pt-wise < 0.05'})
    set(gca,'FontSize',12,'FontName','Arial')
end

end