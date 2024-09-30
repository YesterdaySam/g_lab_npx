function [fhandle] = plot_regional_ADFR(root,varargin)
%% Plots the amplitude by depth of each unit
% Separates colormaps by unit classification (noise/mua/good)
% dot size indicates firing rate (normalized)
% dot color differentiates individual units within colormap/sorting group
%
% Inputs:
%   root        = root object. Must have root.tssync and root.tsb fields
%   varargin
%       unitinds
%       colormap
% Outputs:
%   fhandle = handle to figure
%
% Created 8/21/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

for i = 1:(nargin-1)/2
    nUnits(i) = length(varargin{i});
end

normfr      = root.info.fr / max(root.info.fr) * 500;

fhandle = figure; hold on
set(gcf,'units','normalized','position',[0.4 0.2 0.3 0.6])

for i = 1:(nargin-1)/2
    tmpscatter = scatter(root.info.amp(varargin{i}), root.info.depth(varargin{i}), normfr(varargin{i}), varargin{i+(nargin-1)/2}(1:nUnits(i),:), 'filled');
end

xlim([0 max(root.info.amp)]);
ylim([0 max(root.info.depth)]);
ylabel('Distance from tip (um)'); 
xlabel('Amplitude (uV)')
set(gca,'FontSize',12,'FontName','Arial')

end