function [fhandle] = text2corr(fhandle,textstring,mdlparams,xloc)

arguments
    fhandle
    textstring
    mdlparams
    xloc = 0.4;
end

ylabel(textstring)

ylims = ylim;
xlims = xlim;

text(xlims(2) - xloc*diff(xlims), ylims(2)-.1*diff(ylims), ['R = ' num2str(mdlparams.r, 3)], 'FontSize', 12)
text(xlims(2) - xloc*diff(xlims), ylims(2)-.15*diff(ylims), ['p = ' num2str(mdlparams.p, 3)], 'FontSize', 12)
% text(xlims(2) - xloc*diff(xlims), ylims(2)-.2*diff(ylims), ['slope = ' num2str(mdlparams.b, 3)], 'FontSize', 12)
% try
%     text(xlims(2) - xloc*diff(xlims), ylims(2)-.25*diff(ylims), ['y-int = ' num2str(mdlparams.yint, 3)], 'FontSize', 12)
% catch
% end

set(gca,'FontSize',12,'FontName','Arial')

end