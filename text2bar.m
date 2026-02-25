function [fhandle] = text2bar(fhandle,textstring,tmpP,xloc,npoints)

arguments
    fhandle
    textstring
    tmpP
    xloc = 0.6;
    npoints = 0;
end

figure(fhandle)     % Make fhandle active

ylabel(textstring)
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

text(xlims(2) - xloc*diff(xlims), ylims(2)-.1*diff(ylims), sigStr, 'FontSize', 12)
text(xlims(2) - xloc*diff(xlims), ylims(2)-.15*diff(ylims), pstr, 'FontSize', 12)
if npoints ~= 0
    text(xlims(2) - xloc*diff(xlims), ylims(2)-.2*diff(ylims), ['n = ' npoints], 'FontSize', 12)
end
end