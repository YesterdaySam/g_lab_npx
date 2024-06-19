
binName = 'KW001_05302024_rec_D7_g0_t0.imec0.lf.bin';
binPath = 'D:\Kelton\probe_data\KW001\KW001_05302024_rec_D7_g0\KW001_05302024_rec_D7_g0_imec0\';
meta = SGLX_readMeta.ReadMeta(binName, binPath);

nSamp = floor(str2double(meta.fileTimeSecs) * SGLX_readMeta.SampRate(meta));
dataArray = SGLX_readMeta.ReadBin(0, nSamp, meta, binName, binPath);

%%
dw = 1;
dLineList = 6;

digArray = SGLX_readMeta.ExtractDigital(dataArray, meta, dw, dLineList);
% for i = 1:numel(dLineList)
%     plot(digArray(i,:));
%     hold on
% end

save('KW001_05302024_sync','digArray')