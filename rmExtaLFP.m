function [root] = rmExtaLFP(root,band,spath)

arguments
    root
    band = 2
    spath = pwd
end

root.lfp = root.lfp(root.uPSDMax(band,:),:);

saveRoot(root,spath)

end