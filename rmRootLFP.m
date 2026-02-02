function [root, lfp] = rmRootLFP(root)
%% Removes extraneous LFP channels from root and stores them in lfp var.
%
% Inputs:
%   root = root object with MxNxP root.lfp (channels x time points x shanks) 
%       Must have root.uPSDMax field
% Outputs:
%   root = root object with 1xNxP root.lfp and slimmed lfp variables
%   lfp = lfp object with all lfp data saved
%
% Created 1/15/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

lfp.name     = root.name;
lfp.lfp      = root.lfp;
lfp.fs_lfp   = root.fs_lfp;
lfp.lfpinfo  = root.lfpinfo;
lfp.lfp_tsb  = root.lfp_tsb;
lfp.bands    = root.bands;
lfp.uPSDMax  = root.uPSDMax;
lfp.uPSD     = root.uPSD;

nshanks = numel(unique(root.lfpinfo.lfpShank));

root.lfp     = root.lfp(root.uPSDMax(2,:),:);
root.uPSD    = root.uPSD(:,root.uPSDMax(2,:));
root.lfpinfo = root.lfpinfo(root.uPSDMax(2,:),:);
root.uPSDMax = [repmat(1:nshanks,size(root.bands,2),1)]; % Set uPSD to match updated LFP info

end