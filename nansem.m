function [sem] = nansem(dat)

ndat = length(dat);
tmpstd = std(dat,'omitnan');
sem = tmpstd ./ sqrt(ndat);
end