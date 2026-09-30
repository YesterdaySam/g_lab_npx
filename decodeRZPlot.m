function [errLocs, errMeans] = decodeRZPlot(useUnits,sname,rootFrst,rootLast,sessFrst,sessLast,datF,datN)
% Wrapper to automate decoding FxF, NxN and NxF

tau = 0.5; 
usePyrs = rootFrst.good(useUnits);

for i = 1:sessFrst.nlaps % Leave 1 out X-val
    useIndF = sessFrst.lapInclude & sessFrst.runInds;
    useIndF(sessFrst.lapstt(i):sessFrst.lapend(i)) = false;
    tmpFmap = get_rmap(rootFrst,sessFrst,useUnits,useIndF);

    dcIndF = false(size(useIndF));
    dcIndF(sessFrst.lapstt(i):sessFrst.lapend(i)) = true;
    dcIndF = dcIndF & sessFrst.runInds;

    [tmpfxfDecode(i)] = decodePosBayes(rootFrst,sessFrst,tmpFmap,usePyrs,dcIndF,tau);
end

fxfDecode.rPos = vertcat(tmpfxfDecode.rPos);
fxfDecode.dPos = vertcat(tmpfxfDecode.dPos);
fxfDecode.dMat = horzcat(tmpfxfDecode.dMat);
fxfDecode.dErr = vertcat(tmpfxfDecode.dErr);
fxfDecode.newT = vertcat(tmpfxfDecode.newT);

% [fxfDecode] = decodePosBayes(rootFrst,sessFrst,datF.posfr(useUnits,:),usePyrs,tau);
% [nxnDecode] = decodePosBayes(rootLast,sessLast,datN.posfr(useUnits,:),usePyrs,tau);
% [nxfDecode] = decodePosBayes(rootLast,sessLast,datF.posfr(useUnits,:),usePyrs,tau);

[errLocs, errMeans, fhandle1] = plot_bayesErrHisto(0:0.05:0.95,fxfDecode,nxnDecode,nxfDecode);
fhandle2 = plot_postProbXTime(fxfDecode,[sessFrst.ts(sessFrst.lapstt(10)) sessFrst.ts(sessFrst.lapend(15))]);
fhandle3 = plot_postProbXTime(nxnDecode,[sessLast.ts(sessLast.lapstt(10)) sessLast.ts(sessLast.lapend(15))]);
fhandle4 = plot_postProbXTime(nxfDecode,[sessLast.ts(sessLast.lapstt(10)) sessLast.ts(sessLast.lapend(15))]);

fsave(fhandle1,[sname '_bayesErrHisto'],1,0,0)
fsave(fhandle2,[sname '_bayesPosEst_FxF10-15'],1,0,0)
fsave(fhandle3,[sname '_bayesPosEst_NxN10-15'],1,0,0)
fsave(fhandle4,[sname '_bayesPosEst_NxF10-15'],1,0,0)

save(sname,'fxfDecode','nxnDecode','nxfDecode')
end