function [fhandle] = text2bar(fhandle,labelstr,tmpP,xloc,yloc,col,npoints)

arguments
    fhandle
    labelstr
    tmpP
    xloc = 0.6
    yloc = 0.1
    col = 'k'
    npoints = 0
end

figure(fhandle)     % Make fhandle active

ylabel(labelstr)
ylims = ylim;
xlims = xlim;

pstr = ['p = ' num2str(tmpP, 3)];

if tmpP < 0.001
    sigStr = '***';
    pstr = 'p < 0.001';
elseif tmpP < 0.01
    sigStr = '**'; 
elseif tmpP < 0.05 
    sigStr = '*';
else
    sigStr = 'n.s.';
end

text(xlims(2) - xloc*diff(xlims), ylims(2)-yloc*diff(ylims), sigStr, 'FontSize', 12, 'Color',col)
text(xlims(2) - xloc*diff(xlims), ylims(2)-(yloc + 0.05)*diff(ylims), pstr, 'FontSize', 12, 'Color',col)
if npoints ~= 0
    text(xlims(2) - xloc*diff(xlims), ylims(2)-(yloc + 0.1)*diff(ylims), ['n = ' npoints], 'FontSize', 12)
end
end