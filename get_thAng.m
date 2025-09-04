function [thAngs, thMRLs, thPs, thZs] = get_thAng(thStats)
%% Utility to quickly scrape theta angles from output of plot_thetaMod

for i = 1:max(size(thStats))
    thMRLs(i) = thStats(i).mrl;
    thAngs(i) = thStats(i).ang;
    thPs(i) = thStats(i).p;
    thZs(i) = thStats(i).z;
end

thAngs = rad2deg(thAngs);

end