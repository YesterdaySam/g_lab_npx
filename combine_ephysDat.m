function [] = combine_ephysDat(datT,fname,sdir)
%% 
% Inputs:
%   datT = a table organizing sessions by mouse and recording day
%   fname = save name of the file containing data variables
%   sdir = location to save combined data variables
% Outputs:
%   None (variables saved in-function)
%
% Updated 5/18/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

celID = [];     % [mouseID, recDay, unit ID, layer ID]
useCC = [];     % Outcome of useUnits (in-layer, >0.1Hz, putative Pyr)
lcDat = [];     % [sigSI, trueSI, truePk, trueLc]
lcMap = [];

ct    = 1;

disp('Combining ephys data')

for i = 1:height(datT)

    % === Load data ===
    cd(datT.fpath{i})

    rootfile = dir("*_root.mat");
    load(rootfile.name)
    % sessfile = dir("*_session.mat");
    % load(sessfile.name)
    epochfile = dir("*_dat.mat");
    load(epochfile.name)
    disp(root.name)

    % if ~isfield(sess,'rwdTrials')
    %     sess = importBhvr(datT.fpath{i});
    % end

    % pkExclude = datStruc.truePk < 1;
    useUnits = root.info.fr(root.goodind) > 0.1 & root.info.uType(root.goodind);
    nCCs = length(root.good);

    % === Concatenate recording data ===
    root = get_lyrXtheta(root,datStruc.thetastats);

    useCC = logical([useCC; useUnits]);
    celID = [celID; str2num(datT.mouse{i}(end-2:end))*ones(nCCs,1), ...
        datT.session(i)*ones(nCCs,1), root.good, root.info.lyrID(root.goodind)];

    % === Concatenate SI and Peak data ===
    if ~isfield(datStruc,'sigSI')
        datStruc.sigSI = sum(datStruc.shufSI > datStruc.trueSI,2) / size(datStruc.shufSI,2);
        datStruc = rmfield(datStruc,'sig');
        save([root.name '_ec_pilot_dat'],'datStruc')
    end

    lcDat = [lcDat; datStruc.sigSI, datStruc.trueSI, datStruc.truePk, datStruc.trueLc];
    lcMap = [lcMap; datStruc.posfr];

    ct = ct + 1;

end

cd(sdir)

save(fname,'celID','useCC','lcDat','lcMap')

end