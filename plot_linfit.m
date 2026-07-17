function [fhandle] = plot_linfit(xdat,ydat)

tmpmdl = get_linfit(xdat,ydat);

figure; hold on;
plot(xdat,ydat,'k.')
plot(xdat,tmpmdl.ypred,'r','LineWidth',2)

end