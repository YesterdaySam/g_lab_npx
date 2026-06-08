function [] = combine_bhvrDat(datT,fname,sdir,sessType)
%% 
% Inputs:
%   datT = a table organizing sessions by mouse and recording day
%   fname = save name of the file containing data variables
%   sdir = location to save combined data variables
%   sessType = Regular = 1; RZ Shift = 2; RZ Rand = 3
% Outputs:
%   None (variables saved in-function)
%
% Updated 5/11/26 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

recID = [];     % [mouseID, recDay]
bvDat = [];     % [frstHalf.lckDI, frstHalf.preRZV, lastHalf.lckDI, lastHalf.preRZV

ct    = 1;
% rzPos = [0.9 1.8];

for i = 1:height(datT)

    % === Load data ===
    cd(datT.fpath{i})

    sessfile = dir("*_session.mat");
    load(sessfile.name)
    disp(sess.name)

    if ~isfield(sess,'rwdTrials')
        sess = importBhvr(datT.fpath{i});
    end

    % === Concatenate recording data ===
    recID = [recID; str2num(datT.mouse{i}(end-2:end)), datT.session(i)];

    bvDat(ct).rzPos = [datT.rzloc1(i) datT.rzloc2(i)];

    if sessType == 1
        % === Concatenate Behavior data ===

        [binedges,~,bvDat(ct).lckRMap] = plot_lickpos(sess,0.03,0);
        [~,bvDat(ct).lckDI] = get_lickDiscrim(sess,bvDat(ct).rzPos*100);
        bvDat(ct).ulckDI = mean(bvDat(ct).lckDI,'omitnan');

        bvDat(ct).spVelMap = plot_lap_velCCorr(sess,0.03,0.002,0);
        bvDat(ct).spVCorr = corr(bvDat(ct).spVelMap', mean(bvDat(ct).spVelMap(end-29:end,:))','rows','complete');

    elseif sessType == 2
        try
            datFile = dir('*_dat.mat');
            load(datFile.name)
        catch
            [sessFrst, sessLast] = splitRec(sess);
        end

        [~,~,bvDat(ct).preLMap] = plot_lickpos(sessFrst,0.03,0);
        [~,~,bvDat(ct).pstLMap] = plot_lickpos(sessLast,0.03,0);

        [~,bvDat(ct).preLckDI] = get_lickDiscrim(sessFrst,bvDat(ct).rzPos*100);
        [~,bvDat(ct).pstLckDI] = get_lickDiscrim(sessLast,bvDat(ct).rzPos*100);
        bvDat(ct).uPreLckDI = mean(bvDat(ct).preLckDI,'omitnan');
        bvDat(ct).uPstLckDI = mean(bvDat(ct).pstLckDI,'omitnan');

        bvDat(ct).preNLap = length(sessFrst.valTrials);
        bvDat(ct).pstNLap = length(sessLast.valTrials);
        bvDat(ct).preLapRwd = get_trialSuccess(sessFrst);
        bvDat(ct).pstLapRwd = get_trialSuccess(sessLast);
        bvDat(ct).uPreLapRwd = mean(bvDat(ct).preLapRwd,'omitnan');
        bvDat(ct).uPstLapRwd = mean(bvDat(ct).pstLapRwd,'omitnan');

        bvDat(ct).preSpVelMap = plot_lap_velCCorr(sessFrst,0.03,0.002,0);
        bvDat(ct).vCorPre = corr(bvDat(ct).preSpVelMap', mean(bvDat(ct).preSpVelMap(end-29:end,:))','rows','complete');
        bvDat(ct).pstSpVelMap = plot_lap_velCCorr(sessLast,0.03,0.002,0);
        try
            bvDat(ct).vCorPre = corr(bvDat(ct).pstSpVelMap', mean(bvDat(ct).pstSpVelMap(31:60,:))','rows','complete');
        catch
            bvDat(ct).vCorPre = corr(bvDat(ct).pstSpVelMap', mean(bvDat(ct).pstSpVelMap)','rows','complete');
        end

    elseif sessType == 3

    end

    ct = ct + 1;

end


cd(sdir)

save(fname,'recID','bvDat')

end