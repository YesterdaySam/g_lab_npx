function [fhandle] = plot_group_datXlyr(dat,datInclude,anatDists,mapcol,plotflag)
%% Plots units in datInclude by distance to layer border and layer center
%
% Inputs:
%   dat = Nx1 array of data for each unit in datInclude
%   datInclude = Nx1 binary of included units. Must match size of dat
%   anatDists = Nx2 array of [dist to center, dist to border] for N units
%   mapcol = Nx3 array of colormap values e.g. autumn(N)
%   plotflag = binary
%
% Outputs:
%   fhandle = handle to figure
%
% Created 8/11/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    dat
    datInclude
    anatDists
    mapcol
    plotflag = 1
end

rng(2); % For repeatability in plotting

if plotflag
    fhandle = figure; hold on
    plot([min(anatDists(datInclude,2)) max(anatDists(datInclude,2))], [0 0], 'k--')
    xrand = 0.05*(rand(sum(datInclude),1)-0.5);
    scatter(anatDists(datInclude,2)+xrand, anatDists(datInclude,1),[],dat(datInclude),'filled')
    colormap(mapcol)
    colorbar;
    xlabel('% Distance to CA1/Sub border')
    ylabel('Distance to layer center (um)')
    set(gca,'FontSize',12,'FontName','Arial','YDir','normal')
end

end