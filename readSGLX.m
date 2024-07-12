
binName = 'KW004_06272024_rec_D1_CA1_g0_t0.imec0.lf.bin';
binPath = 'D:\Kelton\probe_data\KW004\KW004_06272024_rec_D1_CA1_g0\KW004_06272024_rec_D1_CA1_g0_imec0\';
meta = SGLX_readMeta.ReadMeta(binName, binPath);

nSamp = floor(str2double(meta.fileTimeSecs) * SGLX_readMeta.SampRate(meta));
dataArray = SGLX_readMeta.ReadBin(0, nSamp, meta, binName, binPath);

%%
% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
% For 3B2 imec data: the sync pulse is stored in line 6.
dLineList = 6;

digArray = SGLX_readMeta.ExtractDigital(dataArray, meta, dw, dLineList);
% for i = 1:numel(dLineList)
%     plot(digArray(i,:));
%     hold on
% end

save('KW004_06272024_sync','digArray')