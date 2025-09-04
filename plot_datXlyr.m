function [fhandle,d2cs,d2bs] = plot_datXlyr(root,dat,datInclude,borderDists,plotflag)
%% Plots units in datInclude by distance to layer border and layer center
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   dat = Nx1 array of data for each unit in datInclude
%   datInclude = Nx1 binary of included units. Must match size of dat
%   borderDists = array of distances to layer border for each shank
%
% Outputs:
%   fhandle = handle to figure
%   d2cs = Nx1 array of distances to layer center
%   d2bs = Nx1 arra of distances to layer border
%
% Created 8/11/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root
    dat
    datInclude
    borderDists
    plotflag = 1
end

d2cs = get_dist2lyrcenter(root);
d2cs = d2cs(root.goodind);

if ~isfield(root.info,'bdrDist')
    root = get_dist2border(root,borderDists);
end

d2bs = root.info.bdrDist(root.goodind);

rng(2); % For repeatability in plotting

if plotflag
    fhandle = figure; hold on
    plot([min(root.info.bdrDist) max(root.info.bdrDist)], [0 0], 'k--')
    xrand = 0.05*(rand(sum(datInclude),1)-0.5);
    scatter(d2bs(datInclude)+xrand, d2cs(datInclude),[],dat(datInclude),'filled')
    colorbar;
    xlabel('% Distance to CA1/Sub border')
    ylabel('Distance to layer center (um)')
    xlim([-1 1])
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end