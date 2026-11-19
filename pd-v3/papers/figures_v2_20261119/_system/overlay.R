# overlay.R - page overlay for the assembler: 8-pt bold panel letters at each block's top-left corner (baseline
# LETTER_BASE mm below the block top) and thin outlines for pending (blank) blocks. Called by assemble.py.
# Usage: Rscript overlay.R <spec.json> <out.pdf>
args<-commandArgs(TRUE);SYS_DIR<-Sys.getenv('FIGV2_SYS');source(file.path(SYS_DIR,'common.R'))
s<-jsonlite::fromJSON(args[1],simplifyVector=FALSE);W<-s$page_w;H<-s$page_h
raw<-tempfile(fileext='.pdf');open_pdf_mm(raw,W,H,bg='transparent')   # transparent: the overlay must not hide the panels
for(bk in s$blanks)grid.draw(gTree(children=gList(blank_grob(bk$w,bk$h)),vp=viewport(x=unit(bk$x,'mm'),y=unit(H-bk$y,'mm'),width=unit(bk$w,'mm'),height=unit(bk$h,'mm'),just=c('left','top'))))
for(l in s$letters)grid.draw(textGrob(l$letter,x=unit(l$x,'mm'),y=unit(H-l$y-LETTER_BASE,'mm'),just=c('left','bottom'),
 gp=gpar(fontfamily=FONT,fontsize=PT_LETTER,fontface='bold',col=INK)))
invisible(dev.off());pdf_crop(raw,args[2],W,H);unlink(raw)
