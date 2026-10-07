function [lmetbl] = get_lmetable(dat,bvgrp,mID,varnames)
% Assumes dat = [condition1; condition2] where length(condition1) = bvgrp(1).n

splitvar = [zeros(length(dat)/2,1); ones(length(dat)/2,1)];
expgrp = zeros(length(dat)/2,1);
expgrp(bvgrp(1).bvInd) = 1;
expgrp = repmat(expgrp,[2,1]);
lmetbl = table(dat,splitvar,expgrp,[mID; mID],'VariableNames',varnames);
end
