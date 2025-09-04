function [ciup,cidn] = get_CI(dat,sigma)

arguments
    dat
    sigma = 1.96;
end

sem = rmmissing(std(dat,'omitnan')/sqrt(size(dat,1)));
ciup = mean(dat,'omitnan') + sem*sigma;
cidn = mean(dat,'omitnan') - sem*sigma;

end