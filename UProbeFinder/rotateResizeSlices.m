function [] = rotateResizeSlices(fpath)

load(fpath)

for slice = 1:length(sSliceData.Slice)
    sSliceData.Slice(slice).RotateAroundDV = 90;
    sSliceData.Slice(slice).ResizeUpDown = 0.25;
    sSliceData.Slice(slice).ResizeLeftRight = 0.25;
end

[tmpDir, tmpF] = fileparts(fpath);
save([fullfile(tmpDir, tmpF), '_trnsfrm.mat'], 'sSliceData')
end