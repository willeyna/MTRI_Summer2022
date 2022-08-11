% pass in a figure to put a cropped pdf version in /forge/nwilley
function getPDF(fig, name)

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'Inches','PaperSize',[pos(3), pos(4)])
print(fig,name,'-dpdf','-r0')

end