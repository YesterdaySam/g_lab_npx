function [depthMap] = probeDepthMap(root)
%%% channel depth map for LFP channel depth alignment

if root.prbType == 'NPX2.0'
    % NPX2.0 assuming use all 4 tip banks
    chanDepth = [0 0 15 15 30 30 45 45 60 60 75 75 90 90 105 105 120 120 135 135 150 150 165 165 180 180 195 195 210 210 225 225 240 240 255 255 270 270 285 285 300 300 315 315 330 330 345 345 ...
                0 0 15 15 30 30 45 45 60 60 75 75 90 90 105 105 120 120 135 135 150 150 165 165 180 180 195 195 210 210 225 225 240 240 255 255 270 270 285 285 300 300 315 315 330 330 345 345 ...
                360 360 375 375 390 390 405 405 420 420 435 435 450 450 465 465 480 480 495 495 510 510 525 525 540 540 555 555 570 570 585 585 600 600 615 615 630 630 645 645 660 660 675 675 690 690 705 705 ...
                360 360 375 375 390 390 405 405 420 420 435 435 450 450 465 465 480 480 495 495 510 510 525 525 540 540 555 555 570 570 585 585 600 600 615 615 630 630 645 645 660 660 675 675 690 690 705 705 ...
                0 0 15 15 30 30 45 45 60 60 75 75 90 90 105 105 120 120 135 135 150 150 165 165 180 180 195 195 210 210 225 225 240 240 255 255 270 270 285 285 300 300 315 315 330 330 345 345 ...
                0 0 15 15 30 30 45 45 60 60 75 75 90 90 105 105 120 120 135 135 150 150 165 165 180 180 195 195 210 210 225 225 240 240 255 255 270 270 285 285 300 300 315 315 330 330 345 345 ...
                360 360 375 375 390 390 405 405 420 420 435 435 450 450 465 465 480 480 495 495 510 510 525 525 540 540 555 555 570 570 585 585 600 600 615 615 630 630 645 645 660 660 675 675 690 690 705 705 ...
                360 360 375 375 390 390 405 405 420 420 435 435 450 450 465 465 480 480 495 495 510 510 525 525 540 540 555 555 570 570 585 585 600 600 615 615 630 630 645 645 660 660 675 675 690 690 705 705];
    bankMap0 = [1:48, 97:144];
    bankMap1 = [49:96, 145:192];
    bankMap2 = [193:240, 289:336];
    bankMap3 = [241:288, 337:384];
    try
        minChan0 = root.info.depth(find(root.info.ch == 0 | root.info.ch == 1, 1)); %Get depth of first channel on bank 0
        if isempty(minChan0); minChan0 = 0; end
    catch
        minChan0 = 0;
    end
    try 
        minChan1 = root.info.depth(find(root.info.ch == 48 | root.info.ch == 49, 1)); %Get depth of first channel on bank 1
        if isempty(minChan1); minChan1 = 0; end
    catch
        minChan1 = 0;
    end
    try 
        minChan2 = root.info.depth(find(root.info.ch == 192 | root.info.ch == 193, 1)); %Get depth of first channel on bank 2
        if isempty(minChan2); minChan2 = 0; end
    catch
        minChan2 = 0;
    end
    try 
        minChan3 = root.info.depth(find(root.info.ch == 240 | root.info.ch == 241, 1)); %Get depth of first channel on bank 3
        if isempty(minChan3); minChan3 = 0; end
    catch
        minChan3 = 0;
    end

    depthMap(bankMap0) = chanDepth(bankMap0) + minChan0;
    depthMap(bankMap1) = chanDepth(bankMap1) + minChan1;
    depthMap(bankMap2) = chanDepth(bankMap2) + minChan2;
    depthMap(bankMap3) = chanDepth(bankMap3) + minChan3;
else
    % NPX1.0 assuming use tip bank
    chanDepth = (0:383/2)*20;
    chanDepth = repmat(chanDepth,[2,1]);
    chanDepth = reshape(chanDepth,[384,1]);
    try
        minChan0 = root.info.depth(find(root.info.ch == 0 | root.info.ch == 1, 1)); %Get depth of first channel on shank 0
    catch
        minChan0 = 0;
    end
    depthMap = chanDepth + minChan0;
end