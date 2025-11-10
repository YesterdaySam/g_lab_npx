function [mdlparams] = get_linfit(xdat,ydat)

if size(xdat,1) > size(xdat,2)
    xdat = xdat';
end
if size(ydat,1) > size(ydat,2)
    ydat = ydat';
end

mdl = fitlm(xdat, ydat);
ypreds = predict(mdl, xdat');
[r,p] = corrcoef(ydat,ypreds);
r = r(2,1);
p = p(2,1);
b = mdl.Coefficients{2,1}; 

mdlparams.r = r;
mdlparams.p = p;
mdlparams.b = b;
mdlparams.yint = predict(mdl,0);
mdlparams.ypred = ypreds;
end