function save_pdf(f,filename)
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,"pdf_figures/" + filename,'-vector','-dpdf')
end

