function [] = fsave(fhandle,sname,pngF,svgF,figF)
% Saves a variety of images per figure 

arguments
    fhandle
    sname
    pngF = 1
    svgF = 1
    figF = 0
end

if pngF
    saveas(fhandle,sname,'png')
end
if svgF
    saveas(fhandle,sname,'svg')
end
if figF
    saveas(fhandle,sname,'fig')
end
end