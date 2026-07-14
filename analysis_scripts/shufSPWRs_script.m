%% Run shuffles just for SPWRs after LFP channels fixed

spaths       = {'D:\Data\Kelton\analyses\KW106\KW106_06252026_rec_D2_RLat1',...
                'D:\Data\Kelton\analyses\ZM035\ZM035_06112026_rec_D1_LLat1',...
                'D:\Data\Kelton\analyses\ZM035\ZM035_06122026_rec_D2_LLat2'};
datpaths     = {'D:\Data\Kelton\probe_data\KW106\KW106_06252026_rec_D2_RLat1_g0',...
                'D:\Data\Kelton\probe_data\ZM035\ZM035_06112026_rec_D1_LLat1_g0',...
                'D:\Data\Kelton\probe_data\ZM035\ZM035_06122026_rec_D2_LLat2_g0'};
splitRec     = [0,1,2];       % 0 = no split; 1 = RZ shift; 2 = RZ Rand
doUnitParams = true;
doShuffles   = true;
nShuf        = 250;

%%

for j = 1:length(spaths)
    spath = spaths{j};
    datpath = datpaths{j};

    strs = split(spath,'\');
    disp(['Processing SPWRs for ' strs{end}])

    %% Load existing files

    cd(spath)
    
    clear datStruc root sess
    sessfile = dir("*_session.mat");
    load(sessfile.name)
    rootfile = dir("*_root.mat");
    load(rootfile.name)
    datfile = dir("*_dat.mat");
    load(datfile.name)

    try
        if doUnitParams
            if splitRec(j) == 0
                datStruc = get_unitParams(root,root.good,sess,datStruc,false,false,false,true,false,0,false);
            end

            if splitRec(j) > 0
                rootFrst.ripRef = root.ripRef;
                rootLast.ripRef = root.ripRef;
            end
            
            if splitRec(j) == 1
                frstHalf = get_unitParams(rootFrst,root.good,sessFrst,frstHalf,false,false,false,true,false,0,false);
                lastHalf = get_unitParams(rootLast,root.good,sessLast,lastHalf,false,false,false,true,false,0,false);
            elseif splitRec(j) == 2
                frstHalf = get_unitParams(rootFrst,root.good,sessFrst,frstHalf,false,false,false,true,false,0,false);
                lastHalf = get_unitParams(rootLast,root.good,sessLast,lastHalf,false,false,false,true,false,0,false);
            end
        end

        disp('Added SPWR response to units')
    catch
        disp('Failed to add SPWR response to units')
    end

    %% Get shuffle unit parameters

    try
        if doShuffles
            if splitRec(j) == 0
                datStruc = get_shufParams(root,root.good,sess,datStruc,nShuf);
            end

            if splitRec(j) == 1
                frstHalf = get_shufParams(rootFrst,root.good,sessFrst,frstHalf,nShuf,false,true);
                lastHalf = get_shufParams(rootLast,root.good,sessLast,lastHalf,nShuf,false,true);
            elseif splitRec(j) == 2
                frstHalf = get_shufParams(rootFrst,root.good,sessFrst,frstHalf,nShuf,false,true);
                lastHalf = get_shufParams(rootLast,root.good,sessLast,lastHalf,nShuf,false,true);
            end
        end

        disp('Ran SPWR shuffles')
    catch
        disp('Failed to do SPWR shuffles')
    end

    %% Save dat strucs

    try
        if splitRec(j) == 0
            save([root.name '_dat'],'datStruc')
        else
            save([root.name '_dat'],'frstHalf','lastHalf','rootFrst','rootLast','sessFrst','sessLast')
        end

        disp('Saved data structure')
    catch
        disp('Failed to save data structure')
    end

end