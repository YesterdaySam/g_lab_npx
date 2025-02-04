function [ccdepth,ccshank,ccregion] = getDepthByID(root,cc)

ind = find(root.info.cluster_id == cc);

ccdepth = root.info.depth(ind);
ccshank = root.info.shankID(ind);

try
    ccregion = root.info.region(ind);
catch
    disp("No region information in root.info")
end

end