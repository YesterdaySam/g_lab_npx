function [ccdepth,ccshank] = getDepthByID(root,cc)

ind = find(root.info.cluster_id == cc);

ccdepth = root.info.depth(ind);
ccshank = root.info.shankID(ind);

end