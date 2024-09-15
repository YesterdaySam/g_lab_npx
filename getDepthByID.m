function [ccdepth] = getDepthByID(root,cc)

ind = find(root.info.cluster_id == cc);

ccdepth = root.info.depth(ind);

end