function [rmap] = get_rmap(root,sess,useUnits,useInds)

arguments
    root
    sess
    useUnits
    useInds = sess.runInds & sess.lapInclude
end

if islogical(useUnits)
    useUnits = root.good(useUnits);
end

nCCs = length(useUnits);

rmap = nan(nCCs,round(185/5));

for i = 1:length(useUnits)
    [~,~,~,~,~,~,rmap(i,:),~] = get_SI(root,useUnits(i),sess,useInds);
end

end