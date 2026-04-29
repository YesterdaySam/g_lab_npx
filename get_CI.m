function [ciup,cidn] = get_CI(dat,sigma)

arguments
    dat
    sigma = 1.96;
end

sem = rmmissing(std(dat,'omitnan')/sqrt(size(dat,1)));
ciup = rmmissing(mean(dat,'omitnan')) + sem*sigma;
cidn = rmmissing(mean(dat,'omitnan')) - sem*sigma;

end