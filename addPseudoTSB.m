function [root] = addPseudoTSB(rootdir, sdir)
%% Adds root.tsb and lfp.tsb for recordings with no behavior
% Not extensivley tested
% Created 4/24/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

% pseudoFS = 2500;    % Hard coded to match behavior

cd(rootdir)
rootfile = dir('*_root*');
rootstruc = load(rootfile.name); root = rootstruc.root; clear rootstruc;

cd(sdir)

%% Remove double spikes < 1ms apart for each cluster

allcl = unique(root.cl);
badInds = zeros(length(root.ts),1);

for i = 1:height(root.info)
    tmpspks = root.ts(root.cl == allcl(i));
    tmpDoubles = find((diff(tmpspks) < 0.001))+1;   % Index of second spike under 1ms
    for j = 1:length(tmpDoubles)
        tmpInd = root.ts == tmpspks(tmpDoubles(j)) & root.cl == allcl(i);
        badInds(tmpInd) = 1;
    end
end

badInds = logical(badInds);
root.ts(badInds) = [];
root.cl(badInds) = [];

%% Naively assign ts to tsb

root.tsb = root.ts;
root.lfp_tsb = (1:size(root.lfp,2))./root.fs_lfp;

%% Save

save([root.name, '_root'],'root','-v7.3')

end
