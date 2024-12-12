function [root] = batch_get_si(root,sess)
%% Gets Spatial Info of each unit and assigns to root.info.spatial_info
%
% Inputs:
%   root = root object. Must have root.tssync and root.tsb fields
%   sess = session object from importBhvr.mat
%
% Outputs:
%   root = root updated with root.info.spatial_info
%
% Created 12/9/24 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

SI = zeros(1,length(root.good));
nUnits = height(root.info);

for i = 1:nUnits
    cc = root.lb{i,1};
    if root.info.n_spikes(i) < 200
        SI(i) = 0;
    else
    [SI(i)] = get_SI(root,cc,sess);
    end
end

root.info.spatial_info = SI';
end